@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-remapping-ingest', version = '0.6.10-SNAPSHOT')

import com.mongodb.MongoCursorNotFoundException
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.boot.CommandLineRunner
import org.springframework.context.annotation.Import
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import org.springframework.stereotype.Component

import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.remapping.ingest.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.remapping.ingest.parameters.InputParameters

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

@Component
@Import(value=[InputParametersConfiguration.class, MongoConfiguration.class])
class FixRemediation implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(FixRemediation.class)

    @Autowired
    private MongoTemplate mongoTemplate

    @Autowired
    private InputParameters inputParameters

    // Can use the properties file from remapping ingestion, so we have these parameters
    String sourceAssembly = inputParameters.getRemappedFrom()
    String targetAssembly = inputParameters.getAssemblyAccession()

    void run(String... args) {
        int numSplitOperations = 0
        int batchIndex = 0
        String lastSeenID = null
        while (true) {
            ImmutablePair<List<? extends SubmittedVariantOperationEntity>, String> splitOperationsAndLastSeenID =
                    getNextBatchOfSplitOperations(lastSeenID)
            if (splitOperationsAndLastSeenID != null) {
                List<? extends SubmittedVariantEntity> sves =
                        getNextBatchOfSplitEvaSVEs(splitOperationsAndLastSeenID.left)
                if (sves != null) {
                    removeHashesFromDbsnp(sves)
                } else {
                    logger.warn("Got split operations but found no EVA SVEs")
                }
                numSplitOperations += splitOperationsAndLastSeenID.left.size()
                lastSeenID = splitOperationsAndLastSeenID.right
                logger.info("Processed " + numSplitOperations + " split operations...")
                batchIndex += 1
            } else {
                break
            }
        }
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    ImmutablePair<List<? extends SubmittedVariantOperationEntity>, String> getNextBatchOfSplitOperations(
            String lastSeenID) {
        // Get split operations in source assembly
        Query queryToGetNextBatchOfSVO = query(where("_id").regex("^SS_SPLIT_.+")
                .and("inactiveObjects.seq").is(sourceAssembly))
        if (Objects.nonNull(lastSeenID)) {
            queryToGetNextBatchOfSVO.addCriteria(where("_id").gt(lastSeenID))
        }
        queryToGetNextBatchOfSVO = queryToGetNextBatchOfSVO.with(Sort.by(Sort.Direction.ASC, "_id")).limit(1000)

        List<? extends SubmittedVariantOperationEntity> result =
                this.mongoTemplate.find(queryToGetNextBatchOfSVO, DbsnpSubmittedVariantOperationEntity.class)
        if (result.size() > 0) {
            return new ImmutablePair<>(result, result.get(result.size() - 1).getId())
        }
        return null
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    List<? extends SubmittedVariantEntity> getNextBatchOfSplitEvaSVEs(
            List<? extends SubmittedVariantOperationEntity> operations) {
        List<Long> accessions = operations.collect{ it.getSplitInto() }
        // Get split EVA SVEs from target assembly
        Query queryToGetNextBatchOfSVE = query(where("accession").in(accessions)
                .and("seq").is(targetAssembly)
                .and("remappedFrom").exists(true))
        List<? extends SubmittedVariantEntity> result = this.mongoTemplate.find(queryToGetNextBatchOfSVE,
                SubmittedVariantEntity.class)
        if (result.size() > 0) {
            return result
        }
        return null
    }

    void removeHashesFromDbsnp(List<? extends SubmittedVariantEntity> sves) {
        List<String> hashes = sves.collect{it.getHashedMessage()}
        // Delete dbSNP SVEs with same hash from target assembly
        List<DbsnpSubmittedVariantEntity> removed = mongoTemplate.findAllAndRemove(
                query(where("_id").in(hashes).and("seq").is(targetAssembly).and("remappedFrom").exists(true)),
                DbsnpSubmittedVariantEntity.class)
        logger.info("Removed " + removed.size() + "dbsnp SVEs")
    }

}
