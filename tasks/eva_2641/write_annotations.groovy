@Grab(group='uk.ac.ebi.eva', module='eva-pipeline', version='2.8-SNAPSHOT')
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ItemWriter
import org.springframework.batch.item.support.CompositeItemWriter
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.data.mongodb.core.MongoOperations

import uk.ac.ebi.eva.commons.models.metadata.AnnotationMetadata
import uk.ac.ebi.eva.commons.models.mongo.entity.Annotation
import uk.ac.ebi.eva.pipeline.configuration.io.writers.AnnotationCompositeWriterConfiguration
import uk.ac.ebi.eva.pipeline.io.writers.AnnotationInVariantMongoWriter
import uk.ac.ebi.eva.pipeline.io.writers.AnnotationMongoWriter

@Component
@Import(value=[AnnotationCompositeWriterConfiguration.class])
class MainApp implements CommandLineRunner {

    private ItemWriter<List<Annotation>> variantAnnotationItemWriter

    private ItemWriter<List<Annotation>> annotationItemWriter

    private ItemWriter<List<Annotation>> annotationWriter

    @Autowired
    private MongoOperations mongoOperations

    private String vepCacheVersion

    private String vepVersion

    private final String metadataCollectionName = "annotationMetadata_2_0"

    private final String variantCollectionName = "variants_2_0"

    private final String annotationCollectionName = "annotations_2_0"

    private static Logger logger = LoggerFactory.getLogger(MainApp.class)

    void run(String... args) {
        // TODO fill from ...?
        vepVersion = "1"
        vepCacheVersion = "1"

        annotationItemWriter = annotationItemWriter(mongoOperations)
        variantAnnotationItemWriter = variantAnnotationItemWriter(mongoOperations)
        annotationWriter = compositeAnnotationItemWriter()

        List<Annotation> annotations = new ArrayList<>()  // TODO need a different reader
        logger.info("Read {} annotations", annotations.size())

        annotationWriter.write(annotations)
        logger.info("Done writing annotations")

        writeAnnotationMetadata()
        logger.info("Done writing metadata")
    }

    // Replicate annotation writers outside of a step
    // See here for original definitions: https://github.com/EBIvariation/eva-pipeline/tree/c797bfed8999cf1a6752011e200f413896fd34c7/src/main/java/uk/ac/ebi/eva/pipeline/configuration/io/writers
    CompositeItemWriter<List<Annotation>> compositeAnnotationItemWriter() {
        CompositeItemWriter<List<Annotation>> writer = new CompositeItemWriter<>()
        writer.setDelegates(Arrays.asList(annotationItemWriter, variantAnnotationItemWriter))
        return writer;
    }

    ItemWriter<List<Annotation>> variantAnnotationItemWriter(MongoOperations mongoOperations) {
        return new AnnotationInVariantMongoWriter(mongoOperations, variantCollectionName, vepVersion, vepCacheVersion)
    }

    ItemWriter<List<Annotation>> annotationItemWriter(MongoOperations mongoOperations) {
        return new AnnotationMongoWriter(mongoOperations, annotationCollectionName)
    }

    // Replicate metadata tasklet outside of a step
    // See here for original definition: https://github.com/EBIvariation/eva-pipeline/blob/4920e502f415a3ee84d2817e9d2a496f99894c3b/src/main/java/uk/ac/ebi/eva/pipeline/jobs/steps/tasklets/AnnotationMetadataTasklet.java
    void writeAnnotationMetadata() {
        AnnotationMetadata annotationMetadata = new AnnotationMetadata(vepVersion, vepCacheVersion)
        writeUnlessAlreadyPresent(annotationMetadata)
    }

    void writeUnlessAlreadyPresent(AnnotationMetadata annotationMetadata) {
//        long count = mongoOperations.count(new Query(Criteria.byExample(annotationMetadata)), AnnotationMetadata.class,
//                metadataCollectionName);
//        if (count == 0) {
//            mongoOperations.save(annotationMetadata, metadataCollectionName)
//        }
    }

}
