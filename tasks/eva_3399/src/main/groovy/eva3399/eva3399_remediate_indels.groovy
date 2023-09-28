package eva3399


import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.commons.batch.io.UnwindingItemStreamReader
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.remapping.source.Application

import static eva3399.Utils.getCustomFastaAndAssemblyReportPaths
import static eva3399.eva3371_detect_unnormalized_indels.writeChangedVariantsToDev
import static eva3399.eva3380_collect_colliding_ss_hashes.collectCollidingSSHashes
import static eva3399.SourceTargetAssemblies.sourceTargetAssemblies

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script
// 1) persists normalized indels generated in EVA-3371 to the PROD environment
// 2) deletes merged SS
// 3) remaps the persisted normalized indels to the relevant target assemblies
// 4) back-propagates RS to source assemblies
// 5) deprecates orphaned assemblies
def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file to use for database connection", required: true)
cli.devPropertiesFile(args: 1, "Development properties file to use for database connection", required: true)
cli.assemblyAccession(args: 1, "Assembly to analyze", required: true)
cli.remappingOutputDir(args: 1, "Full path to the directory where remapped VCF files are stored", required: true)
cli.remappingConfigFile(args: 1, "Common config file for all remapping processes", required: true)
cli.fastaDir(args: 1, "Top-level directories containing FASTAs for all assemblies", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    throw new Exception("Invalid command line options provided!")
}

String sourceAssembly = options.assemblyAccession
String fastaDir = options.fastaDir
def (customFastaFilePath, customAssemblyReportPath) = getCustomFastaAndAssemblyReportPaths(sourceAssembly, fastaDir)
String outputDirForAssembly = "${options.remappingOutputDir}/GCA_000002035.2"
List<String> targetAssemblies = sourceTargetAssemblies.get(sourceAssembly)


def prodExtractEnv = createFromSpringContext(options.prodPropertiesFile, Application.class,
        ["parameters.assemblyAccession": sourceAssembly, "parameters.fasta": customFastaFilePath,
         "parameters.assemblyReportUrl": "file:" + customAssemblyReportPath,
         "parameters.outputFolder"     : outputDirForAssembly])
def devEnv =
        createFromSpringContext(options.devPropertiesFile, Application.class,
                ["parameters.assemblyAccession": sourceAssembly, "parameters.fasta": customFastaFilePath,
                 "parameters.assemblyReportUrl": "file:" + customAssemblyReportPath,
                 "parameters.outputFolder"     : outputDirForAssembly])
def remediateIndels = new RemediateIndels(sourceAssembly, prodEnv, devEnv)
remediateIndels.runRemediation()

targetAssemblies.each {targetAssembly ->
    def remapNormalizedIndels = new RemapNormalizedIndels(sourceAssembly, targetAssembly, prodExtractEnv,
            options.remappingConfigFile, fastaDir, outputDirForAssembly)
    def remappedVariantsVcf = remapNormalizedIndels.remap()

    def devIngestEnv =
            createFromSpringContext(devPropertiesFile, uk.ac.ebi.eva.remapping.ingest.Application.class,
                    ["parameters.assemblyAccession": targetAssembly, "parameters.vcf": remappedVariantsVcf,
                     "parameters.assemblyReportUrl": "file:" + getCustomFastaAndAssemblyReportPaths(targetAssembly, fastaDir)[1],
                     // No need to worry about these since we are not using any. They are only there to satisfy bean creation
                     "parameters.remappedFrom"     : sourceAssembly, "parameters.loadTo": "EVA",
                     "parameters.remappingVersion" : "", "build.version": "someVersion"])
// Write remapped normalized variants from source to EVA3371 DEV database
    writeChangedVariantsToDev(prodExtractEnv, devIngestEnv)

    def changedVariantsVCFReader =
            new UnwindingItemStreamReader<>(devIngestEnv.springApplicationContext.getBean(VcfReader.class))
    collectCollidingSSHashes(prodExtractEnv, devIngestEnv, changedVariantsVCFReader)
}
// TODO
//BackPropagate
//DeprecateOrphanedRS
