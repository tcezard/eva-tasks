package eva3371

import groovy.cli.picocli.CliBuilder
import groovy.yaml.YamlSlurper
import htsjdk.variant.variantcontext.VariantContext
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemProcessor
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.summary.SubmittedVariantSummaryFunction
import uk.ac.ebi.eva.commons.batch.io.UnwindingItemStreamReader
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.commons.core.models.pipeline.Variant
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.remapping.source.Application
import uk.ac.ebi.eva.commons.core.models.VariantType
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.remapping.ingest.Application as RemappingIngestApplication
import uk.ac.ebi.eva.remapping.source.batch.io.VariantContextWriter
import uk.ac.ebi.eva.remapping.source.configuration.BeanNames as RemappingExtractBeanNames

import java.nio.file.Paths

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script detects unnormalized indels in a given assembly, normalizes them and subsequently
// stores them in the DEV environment in the database eva3371_accession_sharded
def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file to use for database connection", required: true)
cli.devPropertiesFile(args: 1, "Development properties file to use for database connection", required: true)
cli.assemblyAccession(args: 1, "Assembly to analyze", required: true)
cli.fastaDir(args: 1, "Top-level directories containing FASTAs for all assemblies", required: true)
cli.remappingYamlConfigFile(args: 1, "Full path to a remapping YAML config file", required: true)
cli.releaseYamlConfigFile(args: 1, "Full path to a release YAML config file", required: true)
cli.outputDir(args: 1, "Full path to the directory to store output VCF files", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    throw new Exception("Invalid command line options provided!")
}
def logger = LoggerFactory.getLogger(this.class)

Object runProcess(String command) {
    logger.info("Running command: ${command}...")
    def process= ["bash", "-c", command].execute()
    def (out, err) = [new StringBuffer(), new StringBuffer()]
    process.waitForProcessOutput(out, err)
    if(!out.toString().trim().equals("")) logger.info(out.toString())
    if (process.exitValue() != 0) {
        logger.error("Command: $command exited with exit code ${process.exitValue()}!!")
        logger.error("See error messages below: \n" + err.toString())
        throw new Exception(err.toString())
    }
    logger.warn(err.toString())
    logger.info("Command: $command completed successfully!")
    return ["out": out.toString(), "err": err.toString()]
}

String[] getCustomFastaAndAssembly(String assembly, String fastaDir, String customAssemblyScript) {
    // Avoid maxDepth to avoid scanning through a labyrinthe of sub-dirs
    def fastaFilePath = runProcess("find ${fastaDir} -maxdepth 3 -iname " +
            "'${assembly}.fa'| head -1").out.trim()
    def assemblyReportPath = runProcess("find ${fastaDir} -maxdepth 3 -iname " +
            "'${assembly}_assembly_report.txt'| head -1").out.trim()
    if (fastaFilePath.equals("")) throw new Exception("Could not find FASTA file for assembly $assembly!!")
    if (assemblyReportPath.equals("")) throw new Exception("Could not find assembly report for assembly $assembly!!")
    def (customFastaFilePath, customAssemblyReportPath) = [fastaFilePath.replace(".fa", "_custom.fa"),
                                                           assemblyReportPath.replace(".txt", "_custom.txt")]
    runProcess("rm -rf \"${customFastaFilePath}\"")
    runProcess("rm -rf \"${customAssemblyReportPath}\"")
    // Need this to appease the URL-based assembly path resolver here: https://github.com/EBIvariation/eva-accession/blob/33af9bab89219c5071ea77c3c8cd3a1c031a78f3/eva-accession-core/src/main/java/uk/ac/ebi/eva/accession/core/batch/io/AssemblyReportReader.java#L66-L66
    customAssemblyReportPath = "file:${customAssemblyReportPath}"
    // See here: https://github.com/EBIvariation/eva-assembly-ingestion/blob/0086ed78ba5f42b4472b28934d958568f7a5c7ec/eva_assembly_ingestion/nextflow/remap_cluster.nf#L91C25-L91C40
    runProcess("${customAssemblyScript} --assembly-accession ${assembly} " +
            "--fasta-file ${fastaFilePath} --report-file ${assemblyReportPath}")
    return [customFastaFilePath, customAssemblyReportPath]
}

String sortAndIndexVcf (String vcfFile, config) {
    def sortedAndIndexedVcf = vcfFile.replace(".vcf", "_sorted.vcf")
    runProcess("${config.vcfSortScript} -f ${vcfFile} ${sortedAndIndexedVcf}")
    runProcess("${config.executable.bgzip} -f ${sortedAndIndexedVcf}")
    sortedAndIndexedVcf += ".gz"
    runProcess("${config.executable.tabix} -f ${sortedAndIndexedVcf}")
    return sortedAndIndexedVcf
}

String writeAllImpactedSSToVcf(String assembly, EVADatabaseEnvironment env, config, String outputDir) {
    def impactedVariantTypes = [VariantType.INS.toString(), VariantType.DEL.toString(),
                                VariantType.INDEL.toString()]
    def allSveCursors = [sveClass, dbsnpSveClass].collect {
        new RetryableBatchingCursor(where("seq").is(assembly), env.mongoTemplate, it)}
    // We are going to mimic EVA extractor logic (just the processor and the writer parts) here: https://github.com/EBIvariation/eva-accession/blob/33af9bab89219c5071ea77c3c8cd3a1c031a78f3/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/steps/ExportSubmittedVariantsStepConfiguration.java#L79-L79
    // to minimize the pain of writing the SVEs above to VCF
    def sveToVariantContextProcessor =
            env.springApplicationContext.getBean(RemappingExtractBeanNames.SUBMITTED_VARIANT_PROCESSOR,
                    ItemProcessor<SubmittedVariantEntity, VariantContext>.class)
    def outputVCFFileName = Paths.get(outputDir + "/" + "${assembly}_before_norm.vcf").toString()
    // See here: https://github.com/EBIvariation/eva-accession/blob/33af9bab89219c5071ea77c3c8cd3a1c031a78f3/eva-remapping-get-source/src/main/java/uk/ac/ebi/eva/remapping/source/configuration/batch/io/VariantContextWriterConfiguration.java#L34-L34
    def variantContextWriter = new VariantContextWriter(Paths.get(outputVCFFileName), assembly)
    variantContextWriter.open(new ExecutionContext())

    allSveCursors.each{it.each {allSves ->
        def impactedSves = allSves.findAll {sve ->
            def type = VariantType.SNV.toString()
            try {
                type = EVAObjectModelUtils.toClusteredVariant(sve).getType().toString()
            }
            catch (IllegalArgumentException exception) {
            }
            return impactedVariantTypes.contains(type)
        }
        variantContextWriter.write(impactedSves.collect{sveToVariantContextProcessor.process(it)}
                .findAll{Objects.nonNull(it)})
    }}
    variantContextWriter.close()
    outputVCFFileName = sortAndIndexVcf(outputVCFFileName, config)
    return Paths.get(outputVCFFileName).toAbsolutePath().toString()
}

String getNormalizedVcf (String preNormalizedVCF, String customFastaFilePath, config) {
    def postNormalizedVCF = preNormalizedVCF.replace("before_norm", "after_norm")
    runProcess("${config.executable.bcftools} norm -c w -f ${customFastaFilePath} ${preNormalizedVCF} " +
            "-o ${postNormalizedVCF} -Oz")
    runProcess("${config.executable.tabix} -f ${postNormalizedVCF}")
    return postNormalizedVCF
}

String getVcfWithChangedVariants(String preNormalizedVCF, String postNormalizedVCF, config) {
    def changedVCF = postNormalizedVCF.replace("after_norm", "changed_after_norm")
    runProcess("${config.executable.bcftools} isec -w1 -C ${postNormalizedVCF} ${preNormalizedVCF} " +
            "-o ${changedVCF} -Oz")
    return changedVCF
}

void writeChangedVariantsToDev(EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv) {
    // Use components from here: https://github.com/EBIvariation/eva-accession/blob/cfaa462332049730dcf4371f5b374f3af568ba42/eva-accession-clustering/src/main/java/uk/ac/ebi/eva/accession/clustering/configuration/batch/steps/ClusteringFromVcfStepConfiguration.java#L46-L46
    def changedVariantsVCFReader = new UnwindingItemStreamReader<>(devEnv.springApplicationContext.getBean(VcfReader.class))
    changedVariantsVCFReader.open(new ExecutionContext())
    def nextVariant = changedVariantsVCFReader.read()
    def batchSize = 1000
    def variantsPostNorm = new ArrayList<Variant>()
    def numProcessedSoFar = 0
    while(Objects.nonNull(nextVariant)) {
        variantsPostNorm.add(nextVariant)
        nextVariant = changedVariantsVCFReader.read()
        if(variantsPostNorm.size() == batchSize || Objects.isNull(nextVariant)) {
            numProcessedSoFar += variantsPostNorm.size()
            def variantsPostNormByHash = variantsPostNorm.groupBy {
                it.sourceEntries[0].getAttribute("SS_HASH")}
            def hashesToLookInProd = variantsPostNormByHash.keySet()
            def updatedDbsnpAndEvaSves =
                    [dbsnpSveClass, sveClass].collect { collectionClass ->
                    prodEnv.mongoTemplate.find(query(where("_id").in(hashesToLookInProd)),
                            collectionClass).collect {oldSVE ->
                        def variantPostNorm = variantsPostNormByHash.get(oldSVE.hashedMessage)[0]
                        def updatedSV = new SubmittedVariant(oldSVE.referenceSequenceAccession,
                                oldSVE.taxonomyAccession, oldSVE.projectAccession, oldSVE.contig,
                                variantPostNorm.start, variantPostNorm.reference, variantPostNorm.alternate,
                                oldSVE.clusteredVariantAccession, oldSVE.isSupportedByEvidence(),
                                oldSVE.isAssemblyMatch(), oldSVE.isAllelesMatch(), oldSVE.isValidated(),
                                oldSVE.createdDate)
                        def updatedHash = new SubmittedVariantSummaryFunction().andThen(new SHA1HashingFunction()).apply(updatedSV)
                        // Crudely tack on the old SS hash to the start of the remappingId since that is the
                        // least harmless String attribute that allows us to encode extra information
                        // without messing up the object model's integrity
                        // We must be careful with this however when using this info in EVA-3380 to update Indels
                        return new SubmittedVariantEntity(oldSVE.accession, updatedHash, updatedSV, oldSVE.version,
                                oldSVE.remappedFrom, oldSVE.remappedDate,
                                oldSVE.hashedMessage + "_" + oldSVE.remappingId)
                    }}
            devEnv.bulkInsertIgnoreDuplicates(updatedDbsnpAndEvaSves[0], dbsnpSveClass)
            devEnv.bulkInsertIgnoreDuplicates(updatedDbsnpAndEvaSves[1], sveClass)
            variantsPostNorm.clear()
            logger.info("${numProcessedSoFar} variants processed so far...")
        }
    }
    changedVariantsVCFReader.close()
}

String impactedAssembly = options.assemblyAccession
String outputDirForAssembly = options.outputDir + "/" + impactedAssembly
def config = new YamlSlurper().parse(new File(options.remappingYamlConfigFile))
logger.info("Using config:" + config)
// Unfortunately, the VCF sort script is in another config file!!
config.getProperties()["map"]["vcfSortScript"] = new YamlSlurper().parse(new File(options.releaseYamlConfigFile))["vcf-sort-script-path"]
def (customFastaFilePath, customAssemblyReportPath) =
getCustomFastaAndAssembly(impactedAssembly, options.fastaDir, config.executable["custom_assembly"])

runProcess("mkdir -p ${outputDirForAssembly}")
def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class,
        ["parameters.assemblyAccession": impactedAssembly, "parameters.fasta": customFastaFilePath,
         "parameters.assemblyReportUrl": customAssemblyReportPath,
         "parameters.outputFolder"     : outputDirForAssembly])

def preNormalizedVCF = writeAllImpactedSSToVcf(impactedAssembly, prodEnv, config, outputDirForAssembly)
def postNormalizedVCF = getNormalizedVcf(preNormalizedVCF, customFastaFilePath, config)
def changedVariantsVCF = getVcfWithChangedVariants(preNormalizedVCF, postNormalizedVCF, config)
def devEnv =
    createFromSpringContext(options.devPropertiesFile, RemappingIngestApplication.class,
            ["parameters.assemblyAccession": impactedAssembly, "parameters.vcf": changedVariantsVCF,
             "parameters.assemblyReportUrl": customAssemblyReportPath,
             // No need to worry about these since we are not using any. They are only there to satisfy bean creation
             "parameters.remappedFrom": "", "parameters.loadTo": "EVA",
             "parameters.remappingVersion": "", "build.version": "someVersion"])
writeChangedVariantsToDev(prodEnv, devEnv)
