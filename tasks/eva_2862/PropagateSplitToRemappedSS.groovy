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

        Variant variant = vcfReader.read()

        while (variant != null) {
            SubmittedVariantEntity sve = variantProcessor.process(variant)

            if (count < chunkSize) {
                chunk.add(sve)
                count++
            } else {
                processChunk(chunk)
                count = 0
                chunk = new ArrayList<>()
            }
            variant = vcfReader.read()
        }

        // process the last chunk, if it exists
        if (!chunk.isEmpty()) {
            processChunk(chunk)
        }
        vcfReader.close()
    }

    void processChunk(List<? extends SubmittedVariantEntity> chunk) {
        List<? extends SubmittedVariantEntity> remappedSVEWithOldSS = getSVEWithSameHash(chunk)
        if (remappedSVEWithOldSS != null) {
            List<SubmittedVariantEntity> remappedSVEWithUpdatedSS = getUpdatedSVE(remappedSVEWithOldSS, sve)
            mongoTemplate.insert(remappedSVEWithUpdatedSS, SubmittedVariantEntity.class)
            mongoTemplate.findAndRemove(remappedSVEWithOldSS, DbsnpSubmittedVariantEntity.class)
        }
}

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    List<? extends SubmittedVariantEntity> getSVEWithSameHash(List<? extends SubmittedVariantEntity> sves) {
        List<String> hashes = sves.collect{ it.getHashedMessage() }
        Query queryToGetNextBatchOfSVE = query(where("_id").in(hashes))
        List<? extends SubmittedVariantEntity> result = this.mongoTemplate.find(queryToGetNextBatchOfSVE,
                DbsnpSubmittedVariantEntity.class)
        if (result.size() > 0) {
            return result
        }
        return null
    }

    List<? extends SubmittedVariantEntity> getUpdatedSVE(List<? extends SubmittedVariantEntity> sveWithOldSS,
                                                         SubmittedVariantEntity sveWithNewSS) {
        Long newSS = sveWithNewSS.getAccession()
        List<SubmittedVariantEntity> updatedSVE = new ArrayList<>()
        for (SubmittedVariantEntity sve : sveWithOldSS) {
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
