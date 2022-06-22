@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-remapping-get-source', version = '0.6.10-SNAPSHOT')

import com.mongodb.MongoCursorNotFoundException
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemProcessor
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.beans.factory.annotation.Qualifier
import org.springframework.boot.CommandLineRunner
import org.springframework.context.annotation.Import
import org.springframework.data.domain.Sort
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.remapping.source.configuration.batch.io.VariantContextWriterConfiguration
import uk.ac.ebi.eva.remapping.source.configuration.batch.processors.SubmittedVariantsProcessorConfiguration
import uk.ac.ebi.eva.remapping.source.configuration.BeanNames
import uk.ac.ebi.eva.remapping.source.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.remapping.source.batch.io.VariantContextWriter
import uk.ac.ebi.eva.remapping.source.parameters.InputParameters
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

@Component
@Import(value=[VariantContextWriterConfiguration.class, SubmittedVariantsProcessorConfiguration.class,
               InputParametersConfiguration.class, MongoConfiguration.class])
class ExportNewlySplitSS implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(ExportNewlySplitSS.class)

    @Autowired
    private MongoTemplate mongoTemplate

    @Autowired
    private InputParameters inputParameters

    @Autowired
    private ItemProcessor<SubmittedVariantEntity, VariantContext> submittedVariantProcessor

    @Autowired
    @Qualifier(BeanNames.EVA_SUBMITTED_VARIANT_WRITER)
    private VariantContextWriter variantContextWriter

    void run(String... args) {
        // Steps:
        // - read from dbSNP SS operations matching SS_SPLIT_FROM_X_TO_Y
        // - for each Y, get the EVA SS entities with accession Y
        // - process and write variant to VCF using existing processor & writer

        int numSplitOperations = 0
        int batchIndex = 0
        String lastSeenID = null
	    variantContextWriter.open(new ExecutionContext())
        while (true) {
            ImmutablePair<List<? extends SubmittedVariantOperationEntity>, String> splitOperationsAndLastSeenID =
                    getNextBatchOfSplitOperations(lastSeenID)
            if (splitOperationsAndLastSeenID != null) {
                List<? extends SubmittedVariantEntity> sves =
                        getNextBatchOfSplitEvaSVEs(splitOperationsAndLastSeenID.left)
                if (sves != null) {
                    variantContextWriter.write(sves.collect{ submittedVariantProcessor.process(it) })
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
	    variantContextWriter.close()
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    ImmutablePair<List<? extends SubmittedVariantOperationEntity>, String> getNextBatchOfSplitOperations(
            String lastSeenID) {
        Query queryToGetNextBatchOfSVO = query(where("_id").regex("^SS_SPLIT_.+")
                .and("inactiveObjects.seq").is(inputParameters.getAssemblyAccession()))
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
        Query queryToGetNextBatchOfSVE = query(where("accession").in(accessions)
                .and("seq").is(inputParameters.getAssemblyAccession())
                .and("remappedFrom").exists(false))
        List<? extends SubmittedVariantEntity> result = this.mongoTemplate.find(queryToGetNextBatchOfSVE,
                SubmittedVariantEntity.class)
        if (result.size() > 0) {
            return result
        }
        return null
    }

}
