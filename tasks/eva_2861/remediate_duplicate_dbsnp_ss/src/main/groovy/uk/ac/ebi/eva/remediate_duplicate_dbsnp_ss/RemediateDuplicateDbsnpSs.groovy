package uk.ac.ebi.eva.remediate_duplicate_dbsnp_ss

@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-clustering', version = '0.6.10-SNAPSHOT')

import com.mongodb.MongoCursorNotFoundException
import org.apache.commons.lang3.tuple.ImmutablePair
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ItemWriter
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.beans.factory.annotation.Qualifier
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
import uk.ac.ebi.eva.accession.clustering.configuration.BeanNames
import uk.ac.ebi.eva.accession.clustering.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.io.SSSplitWriterConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.listeners.ListenersConfiguration
import uk.ac.ebi.eva.accession.clustering.parameters.InputParameters
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.metrics.configuration.MetricConfiguration

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
    private MongoTemplate mongoTemplate

    @Autowired
    private InputParameters inputParameters

    void run(String... args) {
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

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    ImmutablePair<List<? extends SubmittedVariantEntity>, String> getNextBatchOfDbsnpSVEs(String lastSeenID) {
        Query queryToGetNextBatchOfSS =
                query(where("seq").is(inputParameters.getAssemblyAccession())
                        .and("remappedFrom").exists(false))
        if (Objects.nonNull(lastSeenID)) {
            queryToGetNextBatchOfSS.addCriteria(where("_id").gt(lastSeenID))
        }
        queryToGetNextBatchOfSS = queryToGetNextBatchOfSS.with(Sort.by(Sort.Direction.ASC, "_id")).limit(1000)
        List<? extends SubmittedVariantEntity> result = this.mongoTemplate.find(queryToGetNextBatchOfSS,
                DbsnpSubmittedVariantEntity.class)
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
