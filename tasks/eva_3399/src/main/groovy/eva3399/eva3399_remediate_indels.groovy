package eva3399

import com.mongodb.ReadPreference
import com.mongodb.WriteConcern
import eva3399.*
import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.accession.clustering.Application as ClusteringApplication
import uk.ac.ebi.eva.commons.batch.io.UnwindingItemStreamReader
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.remapping.ingest.Application as RemappingIngestApplication
import uk.ac.ebi.eva.remapping.source.Application as RemappingExtractApplication

import static eva3399.Utils.getCustomFastaAndAssemblyReportPaths
import static eva3399.Utils.runProcess
import static eva3399.eva3371_detect_unnormalized_indels.writeChangedVariantsToDev
import static eva3399.eva3380_collect_colliding_ss_hashes.collectCollidingSSHashes
import static eva3399.SourceTargetAssemblies.*

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script
// 1) persists normalized indels generated in EVA-3371 to the PROD environment
// 2) deletes merged SS
// 3) remaps the persisted normalized indels to the relevant target assemblies
// 4) back-propagates RS to source assemblies
// 5) deprecates orphaned RS
def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file to use for database connection", required: true)
cli.devPropertiesFile(args: 1, "Development properties file to use for database connection", required: true)
cli.assemblyAccession(args: 1, "Assembly to remediate", required: true)
cli.taxonomy(args: 1, "Taxonomy to remediate", required: true)
cli.remappingOutputDir(args: 1, "Full path to the directory where remapped VCF files are stored", required: true)
cli.remappingConfigFile(args: 1, "Common config file for all remapping processes", required: true)
cli.fastaDir(args: 1, "Top-level directories containing FASTAs for all assemblies", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    throw new Exception("Invalid command line options provided!")
}

String assemblyToRemediate = options.assemblyAccession
String fastaDir = options.fastaDir
def (customFastaFilePath, customAssemblyReportPath) = getCustomFastaAndAssemblyReportPaths(assemblyToRemediate, fastaDir)
String outputDirForAssembly = "${options.remappingOutputDir}/${assemblyToRemediate}"
String rsReportPath = "${outputDirForAssembly}/rsReportFile"
runProcess("mkdir -p ${outputDirForAssembly}")
List<String> remappedAssemblies = allRemappedAssemblies.getOrDefault(assemblyToRemediate, [])
// Is this assembly the current target assembly for the taxonomy?
boolean isTargetAssembly = taxonomyTargetAssemblyMap[options.taxonomy].equals(assemblyToRemediate)
Set<String> sourceAssemblies = allRemappedAssemblies.findAll{srcAsm, remappedAsms ->
    remappedAsms.contains(assemblyToRemediate)}.keySet()

def prodExtractEnv = createFromSpringContext(options.prodPropertiesFile,
        RemappingExtractApplication.class,
        ["parameters.assemblyAccession": assemblyToRemediate, "parameters.fasta": customFastaFilePath,
         "parameters.assemblyReportUrl": "file:" + customAssemblyReportPath,
         "parameters.outputFolder"     : outputDirForAssembly])
def prodClusteringEnv = createFromSpringContext(options.prodPropertiesFile,
        ClusteringApplication.class, ["parameters.assemblyAccession": assemblyToRemediate,
                                      "parameters.remappedFrom": "", "parameters.vcf": "",
                                      "parameters.projectAccession": "", "parameters.projects": "",
                                      "parameters.rsReportPath": rsReportPath,
                                      "accessioning.instanceId":"instance-6",
                                      "accessioning.submitted.categoryId":"ss",
                                      "accessioning.clustered.categoryId":"rs",
                                      "accessioning.monotonic.ss.blockSize":100000,
                                      "accessioning.monotonic.ss.blockStartValue":5000000000,
                                      "accessioning.monotonic.ss.nextBlockInterval":1000000000,
                                      "accessioning.monotonic.rs.blockSize":100000,
                                      "accessioning.monotonic.rs.blockStartValue":3000000000,
                                      "accessioning.monotonic.rs.nextBlockInterval":1000000000])
def devEnv =
        createFromSpringContext(options.devPropertiesFile, RemappingExtractApplication.class,
                ["parameters.assemblyAccession": assemblyToRemediate, "parameters.fasta": customFastaFilePath,
                 "parameters.assemblyReportUrl": "file:" + customAssemblyReportPath,
                 "parameters.outputFolder"     : outputDirForAssembly])
prodClusteringEnv.mongoTemplate.setWriteConcern(WriteConcern.MAJORITY)
prodClusteringEnv.mongoTemplate.setReadPreference(ReadPreference.primary())
def remediateIndels = new RemediateIndels(assemblyToRemediate, prodClusteringEnv, devEnv, isTargetAssembly)
remediateIndels.runRemediation()

remappedAssemblies.each {remappedAssembly ->
    def remapNormalizedIndels = new RemapNormalizedIndels(assemblyToRemediate, remappedAssembly, prodExtractEnv,
            options.remappingConfigFile, fastaDir, outputDirForAssembly)
    def remappedVariantsVcf = remapNormalizedIndels.remap()

    def devIngestEnv =
            createFromSpringContext(options.devPropertiesFile, RemappingIngestApplication.class,
                    ["parameters.assemblyAccession": remappedAssembly, "parameters.vcf": remappedVariantsVcf,
                     "parameters.assemblyReportUrl": "file:" + getCustomFastaAndAssemblyReportPaths(remappedAssembly, fastaDir)[1],
                     // No need to worry about these since we are not using any. They are only there to satisfy bean creation
                     "parameters.remappedFrom"     : assemblyToRemediate, "parameters.loadTo": "EVA",
                     "parameters.remappingVersion" : "", "build.version": "someVersion"])
    // Write remapped normalized variants from source to EVA3371 DEV database
    writeChangedVariantsToDev(prodExtractEnv, devIngestEnv)

    def changedVariantsVCFReader =
            new UnwindingItemStreamReader<>(devIngestEnv.springApplicationContext.getBean(VcfReader.class))
    collectCollidingSSHashes(prodExtractEnv, devIngestEnv, changedVariantsVCFReader)
}

if(isTargetAssembly) {
    sourceAssemblies.each {sourceAssembly ->
        new BackPropagateUpdatedRS(assemblyToRemediate, sourceAssembly, prodClusteringEnv).backPropagate()
    }
}

new DeprecateOrphanedRS(assemblyToRemediate, prodClusteringEnv).deprecate()

new RemediationQC(assemblyToRemediate, prodClusteringEnv, devEnv, isTargetAssembly).runQC()
