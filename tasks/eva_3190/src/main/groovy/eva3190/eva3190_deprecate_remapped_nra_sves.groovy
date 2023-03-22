package eva3190

import eva3190.RemappingOutputParams
import eva3190.RemediateDuplicateRemappedSS
import eva3190.DeprecateRemappedNraSS
import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory
import org.springframework.batch.test.MetaDataInstanceFactory
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.commons.batch.io.AggregatedVcfLineMapper
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.commons.core.models.Aggregation
import uk.ac.ebi.eva.commons.core.models.pipeline.Variant
import uk.ac.ebi.eva.remapping.ingest.batch.processors.VariantToSubmittedVariantEntityRemappedProcessor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script deprecates variants remapped from variants due to NRA (novel-reference allele) variants
// generated during remapping - see also https://www.ebi.ac.uk/panda/jira/browse/EVA-2919
def cli = new CliBuilder()
cli.propertiesFile(args: 1, "Properties file to use for deprecation", required: true)
cli.remappedAssembly(args: 1, "Remapped assembly", required: true)
cli.remappedVcfDir(args: 1, "Directory with the remapped variants with NRA alleles", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

// Create a connection to the data source with a maximum pool size 1 i.e., one connection
def dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.remappedAssembly, "spring.datasource.hikari.maximumPoolSize": 1])
def logger = LoggerFactory.getLogger(Application.class)
// Construct variant to SVE processor employed by the remapping ingestion job
def variantToSve = new VariantToSubmittedVariantEntityRemappedProcessor(options.remappedAssembly, "dummySrc", "dummyHash")

def getReaderForVCF = {String vcfFilePath ->
    def vcfFile = new File(vcfFilePath)
    File tempVCFFile = File.createTempFile("tempVCF", ".vcf")
    vcfFile.eachLine { lineFromRemappedVCF ->
        def updatedLine = lineFromRemappedVCF
        if (!(lineFromRemappedVCF.trim().equals("")) && !(lineFromRemappedVCF.startsWith("#"))) {
            // Need this dummy data to satisfy AggregatedVcfLineMapper
            updatedLine += ";AF=1,2,3,4,5,6;AC=1,2,3,4,5,6;"
        }
        tempVCFFile.append(updatedLine + "\n")
    }
    AggregatedVcfLineMapper lineMapper = new AggregatedVcfLineMapper("fileId", "studyId", Aggregation.BASIC, null)
    lineMapper.setIncludeIds(true)
    VcfReader vcfReader = new VcfReader(lineMapper, tempVCFFile)
    vcfReader.open(MetaDataInstanceFactory.createStepExecution().getExecutionContext())
    return vcfReader
}

def getNextBatchOfVariants = {VcfReader reader ->
    def batchSize= 3000
    List<Variant> variantsInBatch = new ArrayList<Variant>()
    (0..<batchSize).each {
        List<Variant> variantsRead = reader.read()
        if (Objects.isNull(variantsRead)) return variantsInBatch
        variantsInBatch.addAll(variantsRead)
    }
    return variantsInBatch
}

def processBatch = {variantsInBatch, svesInRemappedAssembly, svesInSourceAssembly ->
    svesInRemappedAssembly = svesInRemappedAssembly.groupBy{it.hashedMessage}
    svesInSourceAssembly = svesInSourceAssembly.groupBy{it.accession}
    // Process every group of lines in the VCF for a given SS ID
    // and return a list of invalid SS remaps or remaps that are no longer valid
    return variantsInBatch.collect {ssID, List<Variant> variants ->
        logger.info("Processing entries for ${ssID}...")
        def remappedSves = variants.collect{variantToSve.process(it).hashedMessage}.collect{
            svesInRemappedAssembly[it]}.flatten().findAll{Objects.nonNull(it)}
        if (remappedSves.isEmpty()) return null
        def remappingAttributes = variants[0].sourceEntries[0].attributes
        def getRSID = {rsIDString -> return Objects.isNull(rsIDString) ? null:
                Long.parseLong(rsIDString.replace("rs", ""))}
        // Construct remapping attributes from the info field
        def correspondingRemappingAttributes =
                Collections.singletonList(new RemappingOutputParams(options.remappedAssembly, ssID,
                        getRSID(remappingAttributes.getOrDefault("RS", null)),
                        remappingAttributes["st"], remappingAttributes["rac"], true))
        // Sometimes, the SS ID in remapped assembly may not line up with that in the source assembly
        // due to splits in source assembly not propagated to the remapped assembly
        // So to look for source SS IDs, use SS IDs from the remapped assembly along with the ssID in the VCF
        // See description in https://www.ebi.ac.uk/panda/jira/browse/EVA-3190
        def sourceSves = ([ssID] + remappedSves.collect{it.accession}).toSet().collect{
            svesInSourceAssembly[it]}.findAll{Objects.nonNull(it)}.flatten()
        // Leverage code written in EVA-2928 verbatim to determine invalid remaps
        // We won't be repeating the same mistake as EVA-2928 because
        // we are not assuming that  that remapped SS ID and the source SS ID are the same (see above)
        return RemediateDuplicateRemappedSS.determineSSWithInvalidRemaps(remappedSves, correspondingRemappingAttributes,
                sourceSves)
    }.flatten().findAll{Objects.nonNull(it)}.unique{"" + "${it.hashedMessage}_${it.accession}"}
}

def vcfReader = getReaderForVCF("${options.remappedVcfDir}/${options.remappedAssembly}.vcf")
def numEntriesScanned = 0
while (true) {
    def variantsInBatch = getNextBatchOfVariants(vcfReader)
    if (variantsInBatch.size() == 0) break
    // Group all the variants read by SS ID
    variantsInBatch = variantsInBatch.groupBy {it.ids[0].replace("ss", "").toLong()}
    def sveHashesToFindInRemappedAssembly = variantsInBatch.values().flatten().collect{
        variantToSve.process(it).hashedMessage}
    def svesInRemappedAssembly = [sveClass, dbsnpSveClass].collect{
        dbEnv.mongoTemplate.find(query(where("_id").in(sveHashesToFindInRemappedAssembly)
                .and("remappedFrom").exists(true)), it)}.flatten()
    // Sometimes, the SS ID in remapped assembly may not line up with that in the source assembly
    // due to splits in source assembly not propagated to the remapped assembly
    // So to look for source SS IDs, use SS IDs from the remapped assembly along with the ones in the VCF
    // See description in https://www.ebi.ac.uk/panda/jira/browse/EVA-3190
    def sourceAssemblyAccessionsToLookup = variantsInBatch.keySet().plus(
            svesInRemappedAssembly.collect{it.accession}.toSet())
    def sourceAssembliesToLookup = svesInRemappedAssembly.collect{it.remappedFrom}.toSet()
    def svesInSourceAssembly = [sveClass, dbsnpSveClass].collect{
        dbEnv.mongoTemplate.find(query(where("seq").in(sourceAssembliesToLookup)
                .and("accession").in(sourceAssemblyAccessionsToLookup)), it)}.flatten()
    def svesToDeprecate = processBatch(variantsInBatch, svesInRemappedAssembly, svesInSourceAssembly)
    logger.info("Deprecating ${svesToDeprecate.size()} variants...")
    DeprecateRemappedNraSS.deprecateSS(dbEnv, svesToDeprecate)
    numEntriesScanned += variantsInBatch.size()
    logger.info("${numEntriesScanned} entries scanned so far...")
}
