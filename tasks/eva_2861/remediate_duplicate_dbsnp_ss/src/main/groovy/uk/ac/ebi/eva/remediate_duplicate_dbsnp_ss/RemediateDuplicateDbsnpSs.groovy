package uk.ac.ebi.eva.accession.clustering.batch.io

@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-clustering', version = '0.6.10-SNAPSHOT')

import com.mongodb.MongoCursorNotFoundException
import org.apache.commons.lang3.tuple.ImmutablePair
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ItemWriter
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.beans.factory.annotation.Qualifier
import org.springframework.beans.factory.annotation.Value
import org.springframework.boot.CommandLineRunner

import org.springframework.boot.SpringApplication
import org.springframework.boot.autoconfigure.SpringBootApplication
import org.springframework.context.ConfigurableApplicationContext
import org.springframework.context.annotation.Import
import org.springframework.data.domain.Sort
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.accession.clustering.batch.io.SSSplitWriter
import uk.ac.ebi.eva.accession.clustering.configuration.BeanNames
import uk.ac.ebi.eva.accession.clustering.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.io.SSSplitWriterConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.listeners.ListenersConfiguration
import uk.ac.ebi.eva.accession.clustering.parameters.InputParameters
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.metrics.configuration.MetricConfiguration

import EVADatabaseEnvironment

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

@Component
@Import(value=[SSSplitWriterConfiguration.class,
        InputParametersConfiguration.class, MetricConfiguration.class, ListenersConfiguration.class])
class RemediateDuplicateDbsnpSs implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(RemediateDuplicateDbsnpSs.class)

    @Autowired
    @Qualifier(BeanNames.SS_SPLIT_WRITER)
    private ItemWriter<SubmittedVariantEntity> ssSplitWriter

    @Autowired
    private MongoTemplate prodMongoTemplate

    private MongoTemplate devMongoTemplate

    @Autowired
    private InputParameters inputParameters

    @Value('${devenv.properties}')
    private String devEnvPropertiesFile

    void run(String... args) {
        // Just in case a previous process had crashed, it is better to process any outstanding split candidates
        // that were left unprocessed during the previous run
        processOutstandingSplitCandidates()
        this.devMongoTemplate = EVADatabaseEnvironment.parseFrom(devEnvPropertiesFile).mongoTemplate
        int numSSScanned = 0
        int batchIndex = 0
        String lastSeenID = null
        while(true) {
            ImmutablePair<List<? extends SubmittedVariantEntity>, String> sveAndLastSeenID =
                    getNextBatchOfDbsnpSVEs(lastSeenID)
            if (sveAndLastSeenID != null) {
                ssSplitWriter.write(sveAndLastSeenID.left)
                numSSScanned += sveAndLastSeenID.left.size()
                lastSeenID = sveAndLastSeenID.right
                logger.info("Processed " + numSSScanned + " duplicate IDs...")
                batchIndex += 1
            } else {
                break
            }
        }
    }

    void processOutstandingSplitCandidates() {
        int numSSScanned = 0
        int batchIndex = 0
        String lastSeenID = null
        while(true) {
            ImmutablePair<List<SubmittedVariantOperationEntity>, String> svoeAndLastSeenID =
                    getNextBatchOfSplitCandidates(lastSeenID)
            if (svoeAndLastSeenID != null) {
                ((SSSplitWriter)ssSplitWriter).processSplitCandidates(svoeAndLastSeenID.left)
                ((SSSplitWriter)ssSplitWriter).removeSplitCandidates(svoeAndLastSeenID.left)
                numSSScanned += svoeAndLastSeenID.left.size()
                lastSeenID = svoeAndLastSeenID.right
                logger.info("Processed " + numSSScanned + " split candidate operations...")
                batchIndex += 1
            } else {
                break
            }
        }
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    ImmutablePair<List<? extends SubmittedVariantEntity>, String> getNextBatchOfDbsnpSVEs(String lastSeenID) {
        Query queryToGetNextBatchOfSS =
                query(where("seq").is(inputParameters.getAssemblyAccession())
                        .and("remappedFrom").exists(false))
        if (Objects.nonNull(lastSeenID)) {
            queryToGetNextBatchOfSS.addCriteria(where("_id").gt(lastSeenID))
        }
        queryToGetNextBatchOfSS = queryToGetNextBatchOfSS.with(Sort.by(Sort.Direction.ASC, "_id")).limit(1000)
        // For the sake of efficiency, get documents from DEV instance rather than the PROD instance
        // because we have only imported a restricted subset of documents that share duplicate IDs in DEV
        // see LoadDuplicateDbsnpSSToTestDB.groovy
        List<? extends SubmittedVariantEntity> result =
                this.devMongoTemplate.find(queryToGetNextBatchOfSS, DbsnpSubmittedVariantEntity.class)
        if (result.size() > 0) {
            return new ImmutablePair<>(result, result.get(result.size() - 1).getId())
        }
        return null
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    ImmutablePair<List<SubmittedVariantOperationEntity>, String> getNextBatchOfSplitCandidates(String lastSeenID) {
        Query queryToGetNextBatchOfSplitCandidates =
                query(where("_id").regex("^SS_SPLIT_CANDIDATES_${inputParameters.getAssemblyAccession()}_.+"))
        if (Objects.nonNull(lastSeenID)) {
            queryToGetNextBatchOfSplitCandidates.addCriteria(where("_id").gt(lastSeenID))
        }
        queryToGetNextBatchOfSplitCandidates = queryToGetNextBatchOfSplitCandidates.with(Sort.by(Sort.Direction.ASC, "_id")).limit(1000)
        List<SubmittedVariantOperationEntity> result =
                this.prodMongoTemplate.find(queryToGetNextBatchOfSplitCandidates, SubmittedVariantOperationEntity.class)
        if (result.size() > 0) {
            return new ImmutablePair<>(result, result.get(result.size() - 1).getId())
        }
        return null
    }

}


@SpringBootApplication
class RemediateDuplicateDbsnpSsApp  {

    static void main(String[] args) {
        ConfigurableApplicationContext context = SpringApplication.run(RemediateDuplicateDbsnpSsApp.class, args)
        System.exit(SpringApplication.exit(context))
    }

}
