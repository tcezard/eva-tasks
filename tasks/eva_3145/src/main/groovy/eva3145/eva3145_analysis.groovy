package eva3145

import eva3145.AsmMatchMetrics
import groovy.cli.picocli.CliBuilder
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.EVAObjectModelUtils
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

// This script analyzes various metrics about variants with assembly match - see EVA-3145

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning", required: true)
cli.devPropertiesFile(args:1, "Development properties file for accessioning", required: true)
cli.assembly(args:1, "Assembly to carry out the analysis for", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

class EVA3145Analysis {
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv
    String assembly
    Logger scriptLogger
    def cvesOfAsmMatchSSColl, cvesOfAsmMatchSSInOtherAsmColl, cvesMergedIntoAsmMatchColl, cvesMergedFromAsmMatchColl,
            cvesSplitFromAsmMatchColl
    def svesWithAsmMatchColl, svesSharingIDInOtherAsmColl, svesSharingIDWithoutAsmMatchColl, svesSplitFromAsmMatchSSColl,
        opsForSvesMergedIntoAsmMatchSSColl
    Map<String, Long> metricsMap = new HashMap<String, Long>()

    EVA3145Analysis(options) {
        this.prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
        this.devEnv = createFromSpringContext(options.devPropertiesFile, GenericApplication.class)
        this.assembly = options.assembly
        this.scriptLogger = LoggerFactory.getLogger(GenericApplication.class)
        (cvesOfAsmMatchSSColl, cvesOfAsmMatchSSInOtherAsmColl, cvesMergedIntoAsmMatchColl, cvesMergedFromAsmMatchColl,
        cvesSplitFromAsmMatchColl,
        svesWithAsmMatchColl, svesSharingIDInOtherAsmColl, svesSharingIDWithoutAsmMatchColl, svesSplitFromAsmMatchSSColl,
        opsForSvesMergedIntoAsmMatchSSColl) =
        ["cvesOfAsmMatchSS", "cvesOfAsmMatchSSInOtherAsm", "cvesMergedIntoAsmMatch", "cvesMergedFromAsmMatch",
         "cvesSplitFromAsmMatch",
         "svesWithAsmMatch", "svesSharingIDInOtherAsm", "svesSharingIDWithoutAsmMatch", "svesSplitFromAsmMatchSS",
         "svesMergedIntoAsmMatchSSOps"]
        def metrics = ["rs_of_asmmatch_ss", "rs_of_asmmatch_ss_in_other_asm",
                       "rs_merged_into_asmmatch_rs", "rs_merged_from_asmmatch_rs", "rs_split_from_asmmatch_rs",
                       "ss_with_asmmatch", "ss_sharing_id_in_other_asm", "ss_sharing_id_without_asmmatch",
                       "ss_split_from_asmmatch", "ss_merged_into_asmmatch"]
        metrics.each{metricsMap.put(it, 0)}
    }

    def runAsmMatchAnalysis() {
        def dbsnpSveAsmMatchCursor = new RetryableBatchingCursor(where("seq").is(assembly)
                .and("asmMatch").is(false), prodEnv.mongoTemplate, dbsnpSveClass)
        def numSVEsScanned = 0

        dbsnpSveAsmMatchCursor.each { List<DbsnpSubmittedVariantEntity> svesWithAsmMatch ->
            metricsMap["ss_with_asmmatch"] = metricsMap["ss_with_asmmatch"] + svesWithAsmMatch.size().toLong()
            devEnv.bulkInsertIgnoreDuplicates(svesWithAsmMatch, sveClass, svesWithAsmMatchColl)
            // RS with constituent SS that have asmMatch false flag
            collectAsmMatchCVEMetrics(svesWithAsmMatch)
            // SS without asmMatch that share accessions with other SS that do
            collectSVEsWithoutAsmMatch(svesWithAsmMatch)
            // SS split from other SS with assembly match
            collectSVEsSplitFromAsmMatchSS(svesWithAsmMatch)
            // SS merged into other SS with map weight
            collectOpsForSVEsMergedIntoAsmMatchSS(svesWithAsmMatch)
            numSVEsScanned += svesWithAsmMatch.size()
            scriptLogger.info("${numSVEsScanned} dbsnpSVE hashes scanned so far in ${assembly}...")
        }
        metricsMap["rs_of_asmmatch_ss"] += devEnv.mongoTemplate.count(query(where("asm").is(assembly)),
                cvesOfAsmMatchSSColl)
        metricsMap["rs_merged_into_asmmatch_rs"] += devEnv.mongoTemplate.count(query(
                where("inactiveObjects.asm").is(assembly)), cvesMergedIntoAsmMatchColl)
        metricsMap["rs_merged_from_asmmatch_rs"] += devEnv.mongoTemplate.count(query(
                where("inactiveObjects.asm").is(assembly)), cvesMergedFromAsmMatchColl)
        metricsMap["rs_split_from_asmmatch_rs"] += devEnv.mongoTemplate.count(query(where("asm").is(assembly)),
                cvesSplitFromAsmMatchColl)
        metricsMap["ss_sharing_id_without_asmmatch"] += devEnv.mongoTemplate.count(query(where("seq").is(assembly)),
                svesSharingIDWithoutAsmMatchColl)
        metricsMap["ss_split_from_asmmatch"] += devEnv.mongoTemplate.count(query(where("seq").is(assembly)),
                svesSplitFromAsmMatchSSColl)
        metricsMap["ss_merged_into_asmmatch"] += devEnv.mongoTemplate.count(query(
                where("inactiveObjects.seq").is(assembly)), opsForSvesMergedIntoAsmMatchSSColl)

        devEnv.mongoTemplate.save(new AsmMatchMetrics(assembly, metricsMap))
    }

    def collectOpsForSVEsMergedIntoAsmMatchSS = { List<DbsnpSubmittedVariantEntity> svesWithAsmMatch ->
        def opsForSvesMergedIntoAsmMatchSS =
                prodEnv.mongoTemplate.find(query(where("eventType").is("MERGED")
                        .orOperator(where("inactiveObjects.asmMatch").exists(false),
                                where("inactiveObjects.asmMatch").is(true))
                        .and("inactiveObjects.seq").is(assembly)
                        .and("mergeInto").in(svesWithAsmMatch.collect{it.accession})),
                        dbsnpSvoeClass)
        devEnv.bulkInsertIgnoreDuplicates(opsForSvesMergedIntoAsmMatchSS, dbsnpSvoeClass,
                opsForSvesMergedIntoAsmMatchSSColl)
    }

    def collectSVEsSplitFromAsmMatchSS = { List<DbsnpSubmittedVariantEntity> svesWithAsmMatch ->
        def sveHashesSplitFromAsmMatchSS =
                prodEnv.mongoTemplate.find(query(where("_id").regex("SS_SPLIT_FROM_.*")
                        .and("inactiveObjects.seq").is(assembly)
                        .and("accession").in(svesWithAsmMatch.collect{it.accession})),
                        dbsnpSvoeClass).collect{it.inactiveObjects[0].hashedMessage}
        def svesSplitFromAsmMatchSS = prodEnv.mongoTemplate.find(query(where("_id")
                .in(sveHashesSplitFromAsmMatchSS)), sveClass)
        devEnv.bulkInsertIgnoreDuplicates(svesSplitFromAsmMatchSS, sveClass, svesSplitFromAsmMatchSSColl)
    }

    def collectSVEsWithoutAsmMatch = { List<DbsnpSubmittedVariantEntity> svesWithAsmMatch ->
        def ssIDsToLookFor = svesWithAsmMatch.collect{it.accession}
        def ssHashesToAvoid = svesWithAsmMatch.collect{it.hashedMessage}
        def svesSharingAccession = prodEnv.mongoTemplate.find(query(where("accession").in(ssIDsToLookFor)
                .and("_id").nin(ssHashesToAvoid)), dbsnpSveClass)
        devEnv.bulkInsertIgnoreDuplicates(svesSharingAccession.findAll{it.isAssemblyMatch() &&
                it.referenceSequenceAccession.equals(assembly)} , dbsnpSveClass, svesSharingIDWithoutAsmMatchColl)
        def svesSharingIDInOtherAsm = svesSharingAccession.findAll{
            !(it.referenceSequenceAccession.equals(assembly))}
        metricsMap["ss_sharing_id_in_other_asm"] += svesSharingIDInOtherAsm.size().toLong()
        devEnv.bulkInsertIgnoreDuplicates(svesSharingIDInOtherAsm , dbsnpSveClass, svesSharingIDInOtherAsmColl)
    }

    def collectCvesMergedIntoAsmMatchCVEs = { targetCves ->
        def targetCveAccessions = targetCves.collect{it.accession}.toSet()
        def targetCvesAccessionHashCombo = targetCves.collect{it.accession + "_" + it.hashedMessage}.toSet()
        def cvesMergedIntoAsmMatchCVEs = [cvoeClass, dbsnpCvoeClass].collect { collectionClass ->
            prodEnv.mongoTemplate.find(query(where("inactiveObjects.asm").is(this.assembly)
                    .and("eventType").is(EventType.MERGED)
                    .and("mergeInto").in(targetCveAccessions)), collectionClass)
        }.flatten().findAll{targetCvesAccessionHashCombo.contains(it.mergedInto + "_" + it.inactiveObjects[0].hashedMessage)}
        devEnv.bulkInsertIgnoreDuplicates(cvesMergedIntoAsmMatchCVEs, dbsnpCvoeClass, cvesMergedIntoAsmMatchColl)
    }

    def collectCvesMergedFromAsmMatchCVEs = { sourceCves ->
        def sourceCveAccessions = sourceCves.collect{it.accession}.toSet()
        def sourceCvesAccessionHashCombo = sourceCves.collect{it.accession + "_" + it.hashedMessage}.toSet()
        def cvesMergedFromAsmMatchCVEs = [cvoeClass, dbsnpCvoeClass].collect { collectionClass ->
            prodEnv.mongoTemplate.find(query(where("inactiveObjects.asm").is(this.assembly)
                    .and("eventType").is(EventType.MERGED)
                    .and("accession").in(sourceCveAccessions)), collectionClass)
        }.flatten().findAll{sourceCvesAccessionHashCombo.contains(it.accession + "_" + it.inactiveObjects[0].hashedMessage)}
        devEnv.bulkInsertIgnoreDuplicates(cvesMergedFromAsmMatchCVEs, dbsnpCvoeClass, cvesMergedFromAsmMatchColl)
    }

    def collectCvesSplitFromAsmMatchCVEs = { asmMatchRSIDs ->
        def cveHashesSplitFromAsmMatchRS =
                prodEnv.mongoTemplate.find(query(where("eventType").is("RS_SPLIT")
                        .and("inactiveObjects.asm").is(assembly)
                        .and("accession").in(asmMatchRSIDs)), dbsnpCvoeClass)
                        .collect{it.inactiveObjects[0].hashedMessage}
        def cvesSplitFromAsmMatchSS = prodEnv.mongoTemplate.find(query(where("_id")
                .in(cveHashesSplitFromAsmMatchRS)), cveClass)
        devEnv.bulkInsertIgnoreDuplicates(cvesSplitFromAsmMatchSS, cveClass, cvesSplitFromAsmMatchColl)
    }

    def collectAsmMatchCVEMetrics = { List<DbsnpSubmittedVariantEntity> svesWithAsmMatch ->
        def svesWithRS = svesWithAsmMatch.findAll{Objects.nonNull(it.clusteredVariantAccession)}
        def svesWithRSGroupedByRSHash = svesWithRS.groupBy{EVAObjectModelUtils.getClusteredVariantHash(it)}
        def rsHashesToLookFor = svesWithRSGroupedByRSHash.keySet()
        def correspondingCVEs = [cveClass, dbsnpCveClass].collect{collectionClass ->
            prodEnv.mongoTemplate.find(query(where("_id").in(rsHashesToLookFor)),
                    collectionClass)}.flatten()
        def rsHashesFound = correspondingCVEs.collect{it.hashedMessage}.toSet()
        (rsHashesToLookFor - rsHashesFound).each{scriptLogger.error("RS hash ${it} corresponding to " +
                "SS hash ${svesWithRSGroupedByRSHash[it][0].hashedMessage} could not be found for assembly ${assembly}!!")}
        def rsIDsFound = correspondingCVEs.collect{it.accession}

        def cvesInOtherAsmByID = [cveClass, dbsnpCveClass].collect{collectionClass ->
            prodEnv.mongoTemplate.find(query(where("accession").in(rsIDsFound).and("asm").ne(assembly)),
                    collectionClass)}.flatten()

        devEnv.bulkInsertIgnoreDuplicates(correspondingCVEs, dbsnpCveClass, cvesOfAsmMatchSSColl)
        metricsMap["rs_of_asmmatch_ss_in_other_asm"] += cvesInOtherAsmByID.size().toLong()
        devEnv.bulkInsertIgnoreDuplicates(cvesInOtherAsmByID, dbsnpCveClass, cvesOfAsmMatchSSInOtherAsmColl)

        collectCvesMergedIntoAsmMatchCVEs(correspondingCVEs)
        collectCvesMergedFromAsmMatchCVEs(correspondingCVEs)
        collectCvesSplitFromAsmMatchCVEs(rsIDsFound)
    }
}

new EVA3145Analysis(options).runAsmMatchAnalysis()
