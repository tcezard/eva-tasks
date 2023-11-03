package eva3399

import groovy.yaml.YamlSlurper
import htsjdk.variant.variantcontext.VariantContext
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemProcessor
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.remapping.source.batch.io.VariantContextWriter
import uk.ac.ebi.eva.remapping.source.configuration.BeanNames as RemappingExtractBeanNames

import static eva3399.Utils.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static eva3399.eva3371_detect_unnormalized_indels.*

import java.nio.file.Paths

class RemapNormalizedIndels {
    String sourceAssembly
    String targetAssembly
    EVADatabaseEnvironment remappingEnv
    Object remappingConfig
    String fastaDir
    String outputDir

    RemapNormalizedIndels() {}

    RemapNormalizedIndels(String sourceAssembly, String targetAssembly, EVADatabaseEnvironment remappingEnv,
                          String remappingConfigFile, String fastaDir, String outputDir) {
        this.sourceAssembly = sourceAssembly
        this.targetAssembly = targetAssembly
        this.remappingEnv = remappingEnv
        this.fastaDir = fastaDir
        this.outputDir = outputDir
        this.remappingConfig = new YamlSlurper().parse(new File(remappingConfigFile))
    }

    def _extractVcfUsingOpClass = { VariantContextWriter variantContextWriter, ssOpClass ->
        def ssClass = ssOpClass.equals(svoeClass) ? sveClass: dbsnpSveClass
        def opsForNormalizedSS =
            new RetryableBatchingCursor(where("_id").regex("^EVA3399_UPD_LOCUS_${this.sourceAssembly}_.*"),
                    this.remappingEnv.mongoTemplate, ssOpClass)
        // We are going to mimic EVA extractor logic (just the processor and the writer parts) here: https://github.com/EBIvariation/eva-accession/blob/33af9bab89219c5071ea77c3c8cd3a1c031a78f3/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/steps/ExportSubmittedVariantsStepConfiguration.java#L79-L79
        // to minimize the pain of writing the SVEs above to VCF
        def sveToVariantContextProcessor =
                this.remappingEnv.springApplicationContext.getBean(
                        RemappingExtractBeanNames.SUBMITTED_VARIANT_PROCESSOR,
                        ItemProcessor<SubmittedVariantEntity, VariantContext>.class)
        opsForNormalizedSS.each{svoes ->
            def ssHashesToFind = svoes.collect{it.id.split("_")[-1]}
            def normalizedSves = this.remappingEnv.mongoTemplate.find(
                    query(where("_id").in(ssHashesToFind).and("remappedFrom").exists(false)), ssClass)
            variantContextWriter.write(normalizedSves.collect{sveToVariantContextProcessor.process(it)}
                    .findAll{Objects.nonNull(it)})
        }
    }

    def _extractVcfFromMongo = {
        def outputVCFFileName = Paths.get(outputDir + "/" + "${this.sourceAssembly}_src_asm_norm.vcf").toString()
        // See here: https://github.com/EBIvariation/eva-accession/blob/33af9bab89219c5071ea77c3c8cd3a1c031a78f3/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/io/VariantContextWriterConfiguration.java#L34-L34
        def variantContextWriter = new VariantContextWriter(Paths.get(outputVCFFileName), this.sourceAssembly)
        variantContextWriter.open(new ExecutionContext())
        // Use this: https://github.com/EBIvariation/eva-accession/blob/5969a350066655098a51e026663a7b90e4b9028f/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/steps/ExportSubmittedVariantsStepConfiguration.java#L51
        [svoeClass, dbsnpSvoeClass].each{_extractVcfUsingOpClass(variantContextWriter, it)}
        variantContextWriter.close()
        outputVCFFileName = sortAndIndexVcf(outputVCFFileName, this.remappingConfig)
        return Paths.get(outputVCFFileName).toAbsolutePath().toString()
    }

    def _setupRemappingBinaries = {
        def binaryPath = "${this.outputDir}/bin".toString()
        // See here: https://github.com/EBIvariation/eva-submission/blob/b99c3659a6a60a0898efe2dd8a8096f9b38e7d19/eva_submission/nextflow/remap_and_cluster.nf#L167
        runProcess("mkdir -p ${binaryPath}".toString())
        def executables = this.remappingConfig.executable
        runProcess([executables.bcftools, executables.samtools, executables.bedtools, executables.minimap2, executables.bgzip,
                    executables.tabix].collect {"ln -s -f ${it} ${binaryPath}/".toString()}.join(" && "))
        return binaryPath
    }

    def remap = {
        def sourceVCFFileName = _extractVcfFromMongo()
        def remappedVCFFileName = sourceVCFFileName.replace(".vcf", "_remapped_to_${this.targetAssembly}.vcf")
        def binaryPath = _setupRemappingBinaries()
        def (srcFasta, srcAsmRpt) = getCustomFastaAndAssemblyReportPaths(this.sourceAssembly, this.fastaDir)
        def (targetFasta, targetAsmRpt) = getCustomFastaAndAssemblyReportPaths(this.targetAssembly, this.fastaDir)
        // See here: https://github.com/EBIvariation/eva-submission/blob/b99c3659a6a60a0898efe2dd8a8096f9b38e7d19/eva_submission/nextflow/remap_and_cluster.nf#L171C5-L171C26
        runProcess(["PATH=${binaryPath}:\$PATH", "source ${this.remappingConfig.executable.python_activate}",
                    "${this.remappingConfig.executable.nextflow} run ${this.remappingConfig.nextflow.remapping} -resume " +
                            "--oldgenome ${srcFasta} --newgenome ${targetFasta} " +
                            "--vcffile ${sourceVCFFileName} --outfile ${remappedVCFFileName}"].join(" && "))
        return remappedVCFFileName
    }
}
