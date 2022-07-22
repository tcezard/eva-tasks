@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')

import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemReader
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.boot.CommandLineRunner
import org.springframework.context.annotation.Import
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.commons.batch.io.AggregatedVcfLineMapper
import uk.ac.ebi.eva.commons.batch.io.UnwindingItemStreamReader
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.commons.core.models.Aggregation
import uk.ac.ebi.eva.commons.core.models.pipeline.Variant

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

@Component
@Import(value=[MongoConfiguration.class])
class QCRemediation implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(QCRemediation.class)

    @Autowired
    private MongoTemplate mongoTemplate

    void run(String... args) {
        String assembly = args[0]
        String unremappedFile = args[1]  // unmapped VCF from remapping pipeline
        String uningestedFile = args[2]  // tsv of SS / hash that wasn't found in dbSNP collection

        Set<Long> unremappedAccessions = getUnremappedAccessions(unremappedFile)
        logger.info("Done loading unremapped accessions")
        Set<Long> uningestedAccessions = getUningestedAccessions(uningestedFile)
        logger.info("Done loading uningested accessions")

        Query getOperationsQuery = query(where("_id").regex("^SS_SPLIT_.+").and("inactiveObjects.seq").is(assembly))
        List<DbsnpSubmittedVariantOperationEntity> operations = this.mongoTemplate.find(getOperationsQuery,
                DbsnpSubmittedVariantOperationEntity.class)
        List<Long> splitAccessions = operations.collect{ it.getSplitInto() }

        int numProblems = 0

        for (Long accession : splitAccessions) {
            Query findRemappedByAccession = query(where("accession").is(accession).and("remappedFrom").exists(true))
            List<DbsnpSubmittedVariantEntity> dbsnpEntities = this.mongoTemplate.find(findRemappedByAccession,
                    DbsnpSubmittedVariantEntity.class)
            List<SubmittedVariantEntity> evaEntities = this.mongoTemplate.find(findRemappedByAccession,
                    SubmittedVariantEntity.class)

            // Should not be in dbSNP
            if (dbsnpEntities.size() > 0) {
                logger.error("Found " + dbsnpEntities.size() + " dbSNP entries with accession " + accession)
                numProblems++
            }
            // Should be in EVA, or unremapped, or uningested
            if (evaEntities.size() == 0 && !(accession in unremappedAccessions) && !(accession in uningestedAccessions)) {
                logger.error("Accession " + accession +
                        " not found in EVA collection, even though it was successfully remapped and ingested")
                numProblems++
            }
        }
        logger.info("Out of " + splitAccessions.size() + " split accessions, found " + numProblems + " with problems.")
    }


    Set<Long> getUnremappedAccessions(String unremappedFile) {
        File vcfFile = new File(unremappedFile)
        AggregatedVcfLineMapper lineMapper = new AggregatedVcfLineMapper("dummyFile", "dummyStudy", Aggregation.BASIC,
                                                                         null)
        lineMapper.setIncludeIds(true)
        lineMapper.setRequireEvidence(false)
        ItemReader<Variant> reader = new UnwindingItemStreamReader(new VcfReader(lineMapper, vcfFile))
        reader.open(new ExecutionContext())
        Set<Long> accessions = new HashSet<>()
        Variant variant = reader.read()
        while (variant != null) {
            accessions.add(Long.parseLong(variant.getMainId().substring(2)))
            variant = reader.read()
        }
        reader.close()
        return accessions
    }

    Set<Long> getUningestedAccessions(String uningestedFile) {
        Set<Long> accessions = new HashSet<>()
        File file = new File(uningestedFile)
        file.eachLine { accessions.add(Long.parseLong(it.split("\t")[0])) }
        return accessions
    }

}
