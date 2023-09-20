package eva3380

import eva3380.ClashingSSHashes
import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory
import org.springframework.batch.item.ExecutionContext
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.summary.SubmittedVariantSummaryFunction
import uk.ac.ebi.eva.commons.batch.io.UnwindingItemStreamReader
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.commons.core.models.pipeline.Variant
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.remapping.ingest.Application as RemappingIngestApplication

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script detects hash collisions in SS among the normalized indels (generated in EVA-3371) + SS in PROD and stores them
// in the DEV environment eva3371_accession_sharded.clashingSSHashes collection
def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file to use for database connection", required: true)
cli.devPropertiesFile(args: 1, "Development properties file to use for database connection", required: true)
cli.assemblyAccession(args: 1, "Assembly to analyze", required: true)
cli.normalizedVcfDir(args: 1, "Full path to the directory where normalized VCF files are stored", required: true)
cli.fastaDir(args: 1, "Top-level directories containing FASTAs for all assemblies", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    throw new Exception("Invalid command line options provided!")
}
def logger = LoggerFactory.getLogger(this.class)

// Copied verbatim from here: https://github.com/EBIvariation/eva-tasks/blob/27e4a1bda59fb30e8baabd2500119ae9ffacf05d/tasks/eva_3371/src/main/groovy/eva3371/eva3371_detect_unnormalized_indels.groovy#L48
def runProcess = {String command ->
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

def getCustomAssemblyReportPath = {String assembly, String fastaDir ->
    def customAssemblyReportPath = runProcess("find ${fastaDir} -maxdepth 3 -iname ".toString() +
            "'${assembly}_assembly_report_custom.txt'| head -1".toString()).out.trim()
    if (customAssemblyReportPath.equals("")) throw new Exception("Could not find assembly report for assembly $assembly!!")
    return customAssemblyReportPath
}

// Copied verbatim from here: https://github.com/EBIvariation/eva-tasks/blob/27e4a1bda59fb30e8baabd2500119ae9ffacf05d/tasks/eva_3371/src/main/groovy/eva3371/eva3371_detect_unnormalized_indels.groovy#L154
def getNormalizedDbsnpAndEvaSves = { ArrayList<Variant> variantsPostNorm,
                                     EVADatabaseEnvironment prodEnv ->
    def variantsPostNormByHash = variantsPostNorm.groupBy {
        it.sourceEntries[0].getAttribute("SS_HASH")
    }
    def hashesToLookInProd = variantsPostNormByHash.keySet()
    def normalizedDbsnpAndEvaSves =
            [dbsnpSveClass, sveClass].collect {collectionClass ->
                prodEnv.mongoTemplate.find(query(where("_id").in(hashesToLookInProd)),
                        collectionClass).collect { oldSVE ->
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
                }
            }
    return normalizedDbsnpAndEvaSves
}

def getClashingSSRecordsFromDb = {List<SubmittedVariantEntity> svesFromVCF,
                                  EVADatabaseEnvironment devEnv, EVADatabaseEnvironment prodEnv ->
    def svesFromVcfGroupedByHash = svesFromVCF.groupBy{it.hashedMessage}
    def clashingSvesInDevAndProd = [devEnv, prodEnv].collect{env ->
        [sveClass, dbsnpSveClass].collect{collectionClass ->
            env.mongoTemplate.find(query(where("_id").in(svesFromVcfGroupedByHash.keySet())), collectionClass)
    }.flatten()}
    def clashingSSAccessionsInProd = clashingSvesInDevAndProd[1].collect{it.accession}.toSet()
    def prodSSAccessionsInNormalizedPopulation =
            [sveClass, dbsnpSveClass].collect{collectionClass ->
                devEnv.mongoTemplate.find(query(where("accession").in(clashingSSAccessionsInProd).and("seq").is(
                        svesFromVCF.first().referenceSequenceAccession)), collectionClass)
            }.flatten().collect{it.accession}.toSet()
    // Don't include clashing SS, that are also lined up for normalization, as a merge candidate
    // See here for illustration: https://docs.google.com/spreadsheets/d/13go87r2jR__iMS5YxDmJ55TAFEwMD0COZOjKnd-0934/edit#rangeid=122505042
    // This way we know for sure that the SS mergeTarget/mergee won't change when another SS normalization is processed
    clashingSvesInDevAndProd[1].removeIf{prodSSAccessionsInNormalizedPopulation.contains(it.accession)}
    def clashingSvesGroupedByHash =
            clashingSvesInDevAndProd.flatten().toUnique{"" + "${it.hashedMessage}_${it.accession}"}.groupBy{it.hashedMessage}
    def hashesFromClashingSSCollection =
            devEnv.mongoTemplate.find(query(where("_id").in(svesFromVcfGroupedByHash.keySet())),
                    ClashingSSHashes.class).groupBy{it.ssHash}
    return clashingSvesGroupedByHash.collect{ssHash, svesWithHashFromDB ->
        def clashingHashesRecord =
                hashesFromClashingSSCollection.getOrDefault(ssHash,
                        Arrays.asList(new ClashingSSHashes(svesWithHashFromDB))).get(0)
        def existingSSIDsForHash = clashingHashesRecord.clashingSS.collect{it.accession}.toSet()
        def clashingSvesFromVCF =
                svesFromVcfGroupedByHash.get(ssHash).findAll{!existingSSIDsForHash.contains(it.accession)}
        clashingHashesRecord.clashingSS.addAll(clashingSvesFromVCF)
        return clashingHashesRecord
    }.findAll{it.clashingSS.size() > 1}
}

// Check EVA3371 output of normalized indels against normalized VCF to collect colliding SS
// Colliding SS -> SS that were left out from ingestion into the normalized SS collection (in eva3371_accession_sharded)
// during EVA3371 due to hash collision
def collectCollidingSSHashes = {EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv, changedVariantsVCFReader ->
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
            def (dbsnpSvesFromVCF, evaSvesFromVCF) = getNormalizedDbsnpAndEvaSves(variantsPostNorm, prodEnv)
            List<ClashingSSHashes> clashingSSRecords =
                    getClashingSSRecordsFromDb(evaSvesFromVCF + dbsnpSvesFromVCF, devEnv, prodEnv)
            clashingSSRecords.each{devEnv.mongoTemplate.save(it)}
            variantsPostNorm.clear()
            logger.info("${numProcessedSoFar} variants processed so far...")
        }
    }
    changedVariantsVCFReader.close()
}

// this is equivalent to if __name__ == '__main__' in Python
if (this.getClass().getName().equals('eva3380.eva3380_collect_colliding_ss_hashes')) {
    String impactedAssembly = options.assemblyAccession
    def changedVariantsVCF = "${options.normalizedVcfDir}/${impactedAssembly}/" +
            "${impactedAssembly}_changed_after_norm_sorted.vcf.gz"
    def assemblyReportUrl = "file:" + getCustomAssemblyReportPath(impactedAssembly, options.fastaDir)
    def prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
    def devEnv =
            createFromSpringContext(options.devPropertiesFile, RemappingIngestApplication.class,
                    ["parameters.assemblyAccession": impactedAssembly, "parameters.vcf": changedVariantsVCF,
                     "parameters.assemblyReportUrl": assemblyReportUrl,
                     // No need to worry about these since we are not using any. They are only there to satisfy bean creation
                     "parameters.remappedFrom"     : "", "parameters.loadTo": "EVA",
                     "parameters.remappingVersion" : "", "build.version": "someVersion"])
    // Use VCF Reader component from here: https://github.com/EBIvariation/eva-accession/blob/cfaa462332049730dcf4371f5b374f3af568ba42/eva-accession-clustering/src/main/java/uk/ac/ebi/eva/accession/clustering/configuration/batch/steps/ClusteringFromVcfStepConfiguration.java#L46-L46
    def changedVariantsVCFReader = new UnwindingItemStreamReader<>(devEnv.springApplicationContext.getBean(VcfReader.class))
    collectCollidingSSHashes(prodEnv, devEnv, changedVariantsVCFReader)
}
