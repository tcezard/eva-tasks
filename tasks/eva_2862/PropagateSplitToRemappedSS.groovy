@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-remapping-ingest', version = '0.6.10-SNAPSHOT')

import com.mongodb.MongoCursorNotFoundException
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemReader
import org.springframework.batch.item.ItemStreamReader
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.boot.CommandLineRunner
import org.springframework.context.annotation.Import
import org.springframework.dao.DuplicateKeyException
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import org.springframework.stereotype.Component

import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.commons.batch.io.UnwindingItemStreamReader
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.commons.core.models.pipeline.Variant
import uk.ac.ebi.eva.remapping.ingest.configuration.batch.io.VcfReaderConfiguration
import uk.ac.ebi.eva.remapping.ingest.configuration.batch.processors.VariantProcessorConfiguration
import uk.ac.ebi.eva.remapping.ingest.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.remapping.ingest.configuration.RemappingMetadataConfiguration
import uk.ac.ebi.eva.remapping.ingest.parameters.InputParameters
import uk.ac.ebi.eva.remapping.ingest.batch.processors.VariantToSubmittedVariantEntityRemappedProcessor

import java.io.BufferedWriter

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

@Component
@Import(value=[VcfReaderConfiguration.class, VariantProcessorConfiguration.class, RemappingMetadataConfiguration.class,
               InputParametersConfiguration.class, MongoConfiguration.class])
class PropagateSplitToRemappedSS implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(PropagateSplitToRemappedSS.class)

    @Autowired
    private MongoTemplate mongoTemplate

    @Autowired
    private InputParameters inputParameters

    @Autowired
    private VcfReader baseVcfReader

    @Autowired
    private VariantToSubmittedVariantEntityRemappedProcessor variantProcessor
    

    void run(String... args) {
        // Steps:
        // - read from remapped vcf
        // - for each variant, compute the hash and get the corresponding document from dbSNP SVE
        // - update the SS while preserving everything else
        // - insert to EVA SVE
        // - delete originals from dbSNP SVE

        int chunkSize = 1000
        int count = 0
        List<SubmittedVariantEntity> chunk = new ArrayList<>()

        ItemReader<Variant> vcfReader = getVcfItemReader(baseVcfReader)
        vcfReader.open(new ExecutionContext())

        // File to write uningested ss to
        BufferedWriter uningestedWriter = new File('uningested_variants.tsv').newWriter()

        Variant variant = vcfReader.read()

        while (variant != null) {
            SubmittedVariantEntity sve = variantProcessor.process(variant)

            if (count < chunkSize) {
                chunk.add(sve)
                count++
            } else {
                processChunk(chunk, uningestedWriter)
                count = 0
                chunk = new ArrayList<>()
            }
            variant = vcfReader.read()
        }

        // process the last chunk, if it exists
        if (!chunk.isEmpty()) {
            processChunk(chunk, uningestedWriter)
        }
        vcfReader.close()
        uningestedWriter.flush()
        uningestedWriter.close()
    }

    void processChunk(List<? extends SubmittedVariantEntity> chunk, BufferedWriter uningestedWriter) {
        logger.info("Processing " + chunk.size() + " variants...")
        Map<Long, ? extends SubmittedVariantEntity> ssToRemappedSVE = getSVEWithSameHash(chunk, uningestedWriter)
        if (ssToRemappedSVE != null) {
            logger.info("Found " + ssToRemappedSVE.size() + " SVE with same hash")
            List<SubmittedVariantEntity> remappedSVEWithUpdatedSS = getUpdatedSVE(ssToRemappedSVE)
            try {
                mongoTemplate.insert(remappedSVEWithUpdatedSS, SubmittedVariantEntity.class)
            } catch (DuplicateKeyException exception) {
                // Duplicates on insertion either means we are rerunning, or (less likely) there is already a matching
                // SVE in the EVA collection - either way we can safely ignore.
                logger.warn(exception.toString())
            }
            mongoTemplate.findAllAndRemove(query(where("accession").in(ssToRemappedSVE.keySet())),
                    DbsnpSubmittedVariantEntity.class)
        }
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    Map<Long, ? extends SubmittedVariantEntity> getSVEWithSameHash(List<? extends SubmittedVariantEntity> sves,
                                                                   BufferedWriter uningestedWriter) {
        /* Returns map from provided SS to the SVE with the same hash */
        Map<String, Long> hashToSS = sves.collectEntries { sve -> [sve.getHashedMessage(), sve.getAccession()] }
        Query queryToGetNextBatchOfSVE = query(where("_id").in(hashToSS.keySet()))
        List<? extends SubmittedVariantEntity> result = this.mongoTemplate.find(queryToGetNextBatchOfSVE,
                DbsnpSubmittedVariantEntity.class)
        Map<String, Long> hashAndSSNotFound = hashToSS.findAll { hashAndSS ->
            !(hashAndSS.value in result.collect{ it.getAccession() }) }
        hashAndSSNotFound.each { hashAndSS -> uningestedWriter.writeLine(hashAndSS.value + "\t" + hashAndSS.key) }
        if (result.size() > 0) {
            return result.collectEntries { sve -> [hashToSS[sve.getHashedMessage()], sve] }
        }
        logger.info("No SVE with same hash for this chunk!")
        return null
    }

    List<? extends SubmittedVariantEntity> getUpdatedSVE(Map<Long, ? extends SubmittedVariantEntity> newSSToSVEWithOldSS) {
        List<SubmittedVariantEntity> updatedSVE = new ArrayList<>()
        for (entry in newSSToSVEWithOldSS.entrySet()) {
            Long newSS = entry.getKey()
            SubmittedVariantEntity sve = entry.getValue()
            SubmittedVariant variant = new SubmittedVariant(
                sve.getReferenceSequenceAccession(),
                sve.getTaxonomyAccession(),
                sve.getProjectAccession(),
                sve.getContig(),
                sve.getStart(),
                sve.getReferenceAllele(),
                sve.getAlternateAllele(),
                sve.getClusteredVariantAccession())  // preserve RSID from the old SVE
            updatedSVE.add(new SubmittedVariantEntity(
                newSS, 
                sve.getHashedMessage(),
                variant,
                sve.getVersion(),
                sve.getRemappedFrom(), 
                sve.getRemappedDate(),
                sve.getRemappingId()))
        }

        return updatedSVE
    }


    // Avoid step-scope bean shenanigans
    ItemStreamReader<Variant> getVcfItemReader(VcfReader vcfReader) {
        return new UnwindingItemStreamReader<>(vcfReader)
    }

}
