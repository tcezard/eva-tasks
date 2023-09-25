package eva3399

import htsjdk.variant.variantcontext.VariantContext
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemProcessor
import uk.ac.ebi.ampt2d.commons.accession.persistence.mongodb.document.EventDocument
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.remapping.source.batch.io.VariantContextWriter
import uk.ac.ebi.eva.remapping.source.configuration.BeanNames as RemappingExtractBeanNames

import static com.jayway.jsonpath.Criteria.where
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

import java.nio.file.Paths

class RemapNormalizedIndels {
    String sourceAssembly
    String targetAssembly
    EVADatabaseEnvironment remappingEnv
    String outputDir

    RemapNormalizedIndels() {}

    RemapNormalizedIndels(String sourceAssembly, String targetAssembly, EVADatabaseEnvironment remappingEnv,
                          String outputDir) {
        this.sourceAssembly = sourceAssembly
        this.targetAssembly = targetAssembly
        this.remappingEnv = remappingEnv
        this.outputDir = outputDir
    }

    def _extractNormalizedIndelsVcf = {VariantContextWriter variantContextWriter, ssOpClass ->
        def ssClass = ssOpClass.equals(svoeClass) ? sveClass: dbsnpSveClass
        def opsForNormalizedSS =
            new RetryableBatchingCursor(where("_id").regex("EVA3399_UPD_LOCUS_.*"),
                    this.remappingEnv.mongoTemplate, ssOpClass)
        // We are going to mimic EVA extractor logic (just the processor and the writer parts) here: https://github.com/EBIvariation/eva-accession/blob/33af9bab89219c5071ea77c3c8cd3a1c031a78f3/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/steps/ExportSubmittedVariantsStepConfiguration.java#L79-L79
        // to minimize the pain of writing the SVEs above to VCF
        def sveToVariantContextProcessor =
                this.remappingEnv.springApplicationContext.getBean(
                        RemappingExtractBeanNames.SUBMITTED_VARIANT_PROCESSOR,
                        ItemProcessor<SubmittedVariantEntity, VariantContext>.class)
        opsForNormalizedSS.each{it.each{List<SubmittedVariantOperationEntity> svoes ->
            def ssHashesToFind = svoes.collect{it.id.split("_")[5]}
            def normalizedSves = this.remappingEnv.mongoTemplate.find(query(where("_id").in(ssHashesToFind)), ssClass)
            variantContextWriter.write(normalizedSves.collect{sveToVariantContextProcessor.process(it)}
                    .findAll{Objects.nonNull(it)})
        }}
    }

    def remap = {
        def outputVCFFileName = Paths.get(outputDir + "/" + "${this.sourceAssembly}_src_asm_norm.vcf").toString()
        // See here: https://github.com/EBIvariation/eva-accession/blob/33af9bab89219c5071ea77c3c8cd3a1c031a78f3/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/io/VariantContextWriterConfiguration.java#L34-L34
        def variantContextWriter = new VariantContextWriter(Paths.get(outputVCFFileName), this.sourceAssembly)
        variantContextWriter.open(new ExecutionContext())
        // Use this: https://github.com/EBIvariation/eva-accession/blob/5969a350066655098a51e026663a7b90e4b9028f/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/steps/ExportSubmittedVariantsStepConfiguration.java#L51
        [svoeClass, dbsnpSvoeClass].each{_extractNormalizedIndelsVcf(variantContextWriter, it)}
        variantContextWriter.close()
        outputVCFFileName = sortAndIndexVcf(outputVCFFileName, config)
        return Paths.get(outputVCFFileName).toAbsolutePath().toString()
    }
}
