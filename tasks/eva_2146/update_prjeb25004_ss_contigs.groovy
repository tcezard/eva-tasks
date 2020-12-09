@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-pipeline', version = '0.4.5-SNAPSHOT')
@Grab('com.xlson.groovycsv:groovycsv:1.3')
import com.mongodb.client.model.DeleteOneModel
import org.bson.Document
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemStreamReader
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.context.ApplicationContext
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Criteria
import org.springframework.data.mongodb.core.query.Query
import uk.ac.ebi.ampt2d.commons.accession.core.models.AccessionWrapper
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.SubmittedVariantAccessioningConfiguration
import uk.ac.ebi.eva.accession.core.contig.ContigMapping
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.service.nonhuman.SubmittedVariantAccessioningService
import uk.ac.ebi.eva.accession.core.service.nonhuman.eva.SubmittedVariantAccessioningDatabaseService
import uk.ac.ebi.eva.accession.core.summary.SubmittedVariantSummaryFunction
import uk.ac.ebi.eva.accession.pipeline.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.accession.pipeline.configuration.batch.io.VcfReaderConfiguration
import uk.ac.ebi.eva.accession.pipeline.configuration.batch.processors.VariantProcessorConfiguration
import uk.ac.ebi.eva.accession.pipeline.parameters.InputParameters
import uk.ac.ebi.eva.commons.batch.io.CoordinatesVcfLineMapper
import uk.ac.ebi.eva.commons.batch.io.UnwindingItemStreamReader
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.commons.core.models.pipeline.Variant

@Component
@Import(value=[SubmittedVariantAccessioningConfiguration.class, MongoConfiguration.class,
        VariantProcessorConfiguration.class, InputParametersConfiguration.class, VcfReaderConfiguration.class])
class MainApp implements CommandLineRunner {

    @Autowired
    private SubmittedVariantAccessioningService svAccessioningService

    @Autowired
    private SubmittedVariantAccessioningDatabaseService svDBService

    @Autowired
    private ContigMapping contigMapping

    @Autowired
    private InputParameters inputParameters

    @Autowired
    private ApplicationContext executionContext

    @Autowired
    private MongoTemplate mongoTemplate

    private List<Long> ssIDBatch = new ArrayList<>()

    private static Logger logger = LoggerFactory.getLogger(MainApp.class);

    void run(String... args) {
        def reader = vcfReader()
        reader.open(new ExecutionContext())
        try {
            int variantCount = 0
            int updateCount = 0
            while(true) {
                Variant variant = reader.read()
                if (variant == null) break
                variantCount += 1
                ssIDBatch.add(Long.parseLong(variant.ids.iterator().next().replace("ss","")))
                if (variantCount % 1000 == 0) {
                    def (variantsToInsert, hashesToDelete) = getVariantsToOverwrite(ssIDBatch)
                    batchWrite(variantsToInsert, hashesToDelete)
                    logger.info("$variantCount records scanned so far...")
                    ssIDBatch.clear()
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace()
        }
        finally {
            reader.close()
        }
    }

    void batchWrite(variantsToInsert, hashesToDelete) {
        if (variantsToInsert.size() > 0 || hashesToDelete.size() > 0) {
            try {
                svDBService.save(variantsToInsert)
                mongoTemplate.getCollection("submittedVariantEntity").bulkWrite(hashesToDelete)
            }
            catch (Exception ex) {
                ex.printStackTrace()
            }
        }
    }

    List getVariantsToOverwrite(List<Long> ssIDs) {
        def variantsToInsert = new ArrayList<>()
        def hashesToDelete = new ArrayList<>()
        mongoTemplate.query(SubmittedVariantEntity.class).matching(new Query(Criteria.where("accession").in(ssIDs))).stream().forEach
                {oldSV ->
            def oldHash = oldSV.hashedMessage
            def study = oldSV.projectAccession
            def contigSynonyms = contigMapping.getContigSynonyms(oldSV.contig)
            if (study.equals("PRJEB25004") && contigSynonyms != null && !(oldSV.contig.equals(contigSynonyms.genBank))) {
                def newContig = contigSynonyms.genBank
                def newSubmittedVariant = new SubmittedVariant(oldSV.model)
                newSubmittedVariant.setContig(newContig)
                def newHash = new SubmittedVariantSummaryFunction().andThen(new SHA1HashingFunction()).apply(newSubmittedVariant)
                //println "New variant is: $newSubmittedVariant"
                def newSVWrapper = new AccessionWrapper<ISubmittedVariant, String, Long>(oldSV.accession, newHash,
                        newSubmittedVariant, oldSV.version)
                variantsToInsert.add(newSVWrapper)
                hashesToDelete.add(new DeleteOneModel<>(new Document("_id", oldHash)))
            }
        }
        return [variantsToInsert, hashesToDelete]
    }

    //Copied from https://github.com/EBIvariation/eva-accession/blob/6b48f8c235e9bc3a9d2d114e700531d57a6c7b20/eva-accession-pipeline/src/main/java/uk/ac/ebi/eva/accession/pipeline/configuration/batch/steps/CheckSubsnpAccessionsStepConfiguration.java#L60-L60
    //This is because a Step Scope bean can't be accessed outside of a Spring Step without resorting to ridiculous machinations!!
    ItemStreamReader<Variant> vcfReader() throws IOException {
        VcfReader vcfReader = new VcfReader(new CoordinatesVcfLineMapper(), new File(inputParameters.getOutputVcf()))
        return new UnwindingItemStreamReader<>(vcfReader)
    }
}
