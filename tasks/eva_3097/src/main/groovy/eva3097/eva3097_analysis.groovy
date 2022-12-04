package eva3097

import groovy.cli.picocli.CliBuilder
import org.slf4j.Logger
import org.springframework.boot.SpringApplication
import org.springframework.data.mongodb.core.query.Query
import org.slf4j.LoggerFactory
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import eva3097.MapWtMetrics

import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

// This script analyzes various metrics about variants with mapping weight - see EVA-3097

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning", required: true)
cli.devPropertiesFile(args:1, "Development properties file for accessioning", required: true)
cli.assembly(args:1, "Assembly to carry out the analysis for", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

new EVA3097Analysis(options).runMapWtAnalysis(options.assembly)

class EVA3097Analysis {
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv
    Logger scriptLogger

    EVA3097Analysis(options) {
        this.prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
        this.devEnv = createFromSpringContext(options.devPropertiesFile, GenericApplication.class)
        this.scriptLogger = LoggerFactory.getLogger(SpringApplication.class)
    }

    def runMapWtAnalysis(String assembly) {
        def metricsMap = [ss_with_mapwt: 0, ss_mapwt_no_rs_mapwt: 0,
                          ss_wo_mapwt_sharing_mapwt_accessions: 0, ss_split_from_ss_with_mapwt: 0,
                          ss_merged_into_ss_with_mapwt: 0,
                          rs_with_mapwt: 0, rs_mapwt_no_ss_mapwt: 0, rs_wo_mapwt_sharing_mapwt_accessions: 0,
                          rs_split_from_rs_with_mapwt: 0, rs_merged_into_rs_with_mapwt: 0]
        metricsMap = ssMapWtAnalysis(assembly, metricsMap)
        metricsMap = rsMapWtAnalysis(assembly, metricsMap)
        devEnv.mongoTemplate.save(new MapWtMetrics(assembly, metricsMap))
    }

    def ssMapWtAnalysis = {String assembly, Map metricsMap ->
        /*
        Four stats collected:
        1) RS without map weight with constituent SS that do
        2) SS without map weight that share accessions with other SS that do
        3) SS split from other SS with map weight
        4) Non map-weighted SS merged into SS with map weight
         */
        def (cvesOfMapWtSSWithoutMapWtColl, svesWithoutMapWtColl, svesSplitFromMapWtSSColl,
             opsForSvesMergedIntoMapWtSSColl) =
        ["cvesOfMapWtSSWithoutMapWt", "svesWithoutMapWt", "svesSplitFromMapWtSS", "svesMergedIntoMapWtSSOps"]
        def dbsnpSveMapWtCursor = new RetryableBatchingCursor(where("seq").is(assembly).and("mapWeight").exists(true),
                prodEnv.mongoTemplate, dbsnpSveClass)
        def numSVEsScanned = 0

        dbsnpSveMapWtCursor.each { List<DbsnpSubmittedVariantEntity> svesWithMapWt ->
            metricsMap["ss_with_mapwt"] += svesWithMapWt.size().toLong()
            // RS without map weight with constituent SS that do
            collectCVEsWithoutMapWt(assembly, svesWithMapWt, cvesOfMapWtSSWithoutMapWtColl)
            // SS without map weight that share accessions with other SS that do
            collectSVEsWithoutMapWt(assembly, svesWithMapWt, svesWithoutMapWtColl)
            // SS split from other SS with map weight
            collectSVEsSplitFromMapWtSS(assembly, svesWithMapWt, svesSplitFromMapWtSSColl)
            // SS merged into other SS with map weight
            collectOpsForSVEsMergedIntoMapWtSS(assembly, svesWithMapWt, opsForSvesMergedIntoMapWtSSColl)
            numSVEsScanned += svesWithMapWt.size()
            scriptLogger.info("${numSVEsScanned} dbsnpSVE hashes scanned so far in ${assembly}...")
        }
        metricsMap["ss_mapwt_no_rs_mapwt"] += devEnv.mongoTemplate.count(query(where("asm").is(assembly)),
                cvesOfMapWtSSWithoutMapWtColl)
        metricsMap["ss_wo_mapwt_sharing_mapwt_accessions"] += devEnv.mongoTemplate.count(query(where("seq").is(assembly)),
                svesWithoutMapWtColl)
        metricsMap["ss_split_from_ss_with_mapwt"] += devEnv.mongoTemplate.count(query(where("seq").is(assembly)),
                svesSplitFromMapWtSSColl)
        metricsMap["ss_merged_into_ss_with_mapwt"] += devEnv.mongoTemplate.count(query(where("inactiveObjects.seq").is(assembly)),
                opsForSvesMergedIntoMapWtSSColl)
        return metricsMap
    }

    def collectOpsForSVEsMergedIntoMapWtSS = { String assembly, List<DbsnpSubmittedVariantEntity> svesWithMapWt,
                                               String opsForSvesMergedIntoMapWtSSColl ->
        def opsForSvesMergedIntoMapWtSS =
                prodEnv.mongoTemplate.find(query(where("eventType").is("MERGED")
                        .and("inactiveObjects.mapWeight").exists(false)
                        .and("inactiveObjects.seq").is(assembly).and("mergeInto")
                        .in(svesWithMapWt.collect{it.accession})), dbsnpSvoeClass)
        devEnv.bulkInsertIgnoreDuplicates(opsForSvesMergedIntoMapWtSS, dbsnpSvoeClass, opsForSvesMergedIntoMapWtSSColl)
    }

    def collectSVEsSplitFromMapWtSS = {String assembly, List<DbsnpSubmittedVariantEntity> svesWithMapWt,
                                       String svesSplitFromMapWtSSColl ->
        def sveHashesSplitFromMapWtSS =
                prodEnv.mongoTemplate.find(query(where("_id").regex("SS_SPLIT_FROM_.*")
                        .and("inactiveObjects.seq").is(assembly).and("accession")
                        .in(svesWithMapWt.collect{it.accession})), dbsnpSvoeClass).collect{it.inactiveObjects[0].hashedMessage}
        def svesSplitFromMapWtSS = prodEnv.mongoTemplate.find(query(where("_id").in(sveHashesSplitFromMapWtSS)),
                sveClass)
        devEnv.bulkInsertIgnoreDuplicates(svesSplitFromMapWtSS, sveClass, svesSplitFromMapWtSSColl)
    }

    def collectSVEsWithoutMapWt = {String assembly, List<DbsnpSubmittedVariantEntity> svesWithMapWt,
                                   String svesWithoutMapWtColl ->
        def svesWithoutMapWtThatShareAccession =
                prodEnv.submittedVariantAccessioningService.getAllActiveByAssemblyAndAccessionIn(assembly,
                        svesWithMapWt.collect { it.accession })
                        .findAll { Objects.isNull(it.data.mapWeight) }
                        .collect { new SubmittedVariantEntity(it.accession, it.hash, it.data, 1) }
        devEnv.bulkInsertIgnoreDuplicates(svesWithoutMapWtThatShareAccession, dbsnpSveClass, svesWithoutMapWtColl)
    }

    def collectCVEsWithoutMapWt = {String assembly, List<DbsnpSubmittedVariantEntity> svesWithMapWt,
                                   String cvesOfMapWtSSWithoutMapWtColl ->
        def correspondingCVEsWithoutMapWt =
                prodEnv.clusteredVariantAccessioningService.getAllActiveByAssemblyAndAccessionIn(assembly,
                        svesWithMapWt.findAll{!Objects.isNull(it.clusteredVariantAccession)}
                                .collect { it.clusteredVariantAccession })
                        .findAll { Objects.isNull(it.data.mapWeight) }
                        .collect { new ClusteredVariantEntity(it.accession, it.hash, it.data, 1) }
        devEnv.bulkInsertIgnoreDuplicates(correspondingCVEsWithoutMapWt, dbsnpCveClass, cvesOfMapWtSSWithoutMapWtColl)
    }

    def rsMapWtAnalysis = {String assembly, Map metricsMap ->
        /*
        Four stats collected:
        1) RS with map weight with constituent SS that do not
        2) RS without map weight that share accessions with other RS that do
        3) RS split from other RS with map weight
        4) Non map-weighted RS merged into RS with map weight
         */
        def (svesWithoutMapWtUnderMapWtRSColl, cvesWithoutMapWtSharingMapWtRSColl, cvesSplitFromMapWtRSColl,
             opsForCvesMergedIntoMapWtRSColl) =
        ["svesWithoutMapWtUnderMapWtRS", "cvesWithoutMapWtSharingMapWtRS", "cvesSplitFromMapWtRS", "cvesMergedIntoMapWtRSOps"]
        def dbsnpCveMapWtCursor = new RetryableBatchingCursor(where("asm").is(assembly).and("mapWeight").exists(true),
                prodEnv.mongoTemplate, dbsnpCveClass)
        def numCVEsScanned = 0
        dbsnpCveMapWtCursor.each {List<DbsnpClusteredVariantEntity> cvesWithMapWt ->
            metricsMap["rs_with_mapwt"] += cvesWithMapWt.size().toLong()
            // RS with map weight with constituent SS that do not
            collectSVEsWithoutMapWtUnderMapWtRS(assembly, cvesWithMapWt, svesWithoutMapWtUnderMapWtRSColl)
            // RS without map weight sharing accessions with other RS that do
            collectCVEsWithoutMapWtSharingMapWtRS(assembly, cvesWithMapWt, cvesWithoutMapWtSharingMapWtRSColl)
            // RS split from other RS with map weight
            collectCVEsSplitFromMapWtRS(assembly, cvesWithMapWt, cvesSplitFromMapWtRSColl)
            // RS merged into other RS with map weight
            collectOpsForCVEsMergedIntoMapWtSS(assembly, cvesWithMapWt, opsForCvesMergedIntoMapWtRSColl)
            numCVEsScanned += cvesWithMapWt.size()
            scriptLogger.info("${numCVEsScanned} dbsnpCVE hashes scanned so far in ${assembly}...")
        }
        metricsMap["rs_mapwt_no_ss_mapwt"] += devEnv.mongoTemplate.count(query(where("seq").is(assembly)),
                svesWithoutMapWtUnderMapWtRSColl)
        metricsMap["rs_wo_mapwt_sharing_mapwt_accessions"] += devEnv.mongoTemplate.count(query(where("asm").is(assembly)),
                cvesWithoutMapWtSharingMapWtRSColl)
        metricsMap["rs_split_from_rs_with_mapwt"] += devEnv.mongoTemplate.count(query(where("asm").is(assembly)),
                cvesSplitFromMapWtRSColl)
        metricsMap["rs_merged_into_rs_with_mapwt"] += devEnv.mongoTemplate.count(query(where("inactiveObjects.asm").is(assembly)),
                opsForCvesMergedIntoMapWtRSColl)
        return metricsMap
    }

    def collectSVEsWithoutMapWtUnderMapWtRS = {String assembly, List<DbsnpClusteredVariantEntity> cvesWithMapWt,
                                               String svesWithoutMapWtUnderMapWtRSColl ->
        def svesWithoutMapWtUnderMapWtRS =
                prodEnv.submittedVariantAccessioningService.getByClusteredVariantAccessionIn(
                        cvesWithMapWt.collect { it.accession })
                        .findAll{it.data.referenceSequenceAccession.equals(assembly) && Objects.isNull(it.data.mapWeight)}
                        .collect { new SubmittedVariantEntity(it.accession, it.hash, it.data, 1) }
        devEnv.bulkInsertIgnoreDuplicates(svesWithoutMapWtUnderMapWtRS, dbsnpSveClass, svesWithoutMapWtUnderMapWtRSColl)
    }

    def collectCVEsWithoutMapWtSharingMapWtRS = {String assembly, List<DbsnpClusteredVariantEntity> cvesWithMapWt,
                                                 String cvesWithoutMapWtSharingMapWtRSColl ->
        def cvesWithoutMapWtSharingMapWtRS =
                prodEnv.clusteredVariantAccessioningService.getAllActiveByAssemblyAndAccessionIn(assembly,
                        cvesWithMapWt.collect { it.accession })
                        .findAll { Objects.isNull(it.data.mapWeight) }
                        .collect { new ClusteredVariantEntity(it.accession, it.hash, it.data, 1) }
        devEnv.bulkInsertIgnoreDuplicates(cvesWithoutMapWtSharingMapWtRS, dbsnpCveClass, cvesWithoutMapWtSharingMapWtRSColl)
    }

    def collectCVEsSplitFromMapWtRS = {String assembly, List<DbsnpClusteredVariantEntity> cvesWithMapWt,
                                       String cvesSplitFromMapWtRSColl ->
        def cveHashesSplitFromMapWtSS =
                prodEnv.mongoTemplate.find(query(where("eventType").is("RS_SPLIT")
                        .and("inactiveObjects.asm").is(assembly).and("accession")
                        .in(cvesWithMapWt.collect{it.accession})), dbsnpCvoeClass).collect{it.inactiveObjects[0].hashedMessage}
        def cvesSplitFromMapWtSS = prodEnv.mongoTemplate.find(query(where("_id").in(cveHashesSplitFromMapWtSS)),
                cveClass)
        devEnv.bulkInsertIgnoreDuplicates(cvesSplitFromMapWtSS, cveClass, cvesSplitFromMapWtRSColl)
    }

    def collectOpsForCVEsMergedIntoMapWtSS = { String assembly, List<DbsnpClusteredVariantEntity> cvesWithMapWt,
                                               String opsForCvesMergedIntoMapWtRSColl ->
        def opsForCvesMergedIntoMapWtRS =
                prodEnv.mongoTemplate.find(query(where("eventType").is("MERGED")
                        .and("inactiveObjects.mapWeight").exists(false)
                        .and("inactiveObjects.asm").is(assembly).and("mergeInto")
                        .in(cvesWithMapWt.collect{it.accession})), dbsnpCvoeClass)
        devEnv.bulkInsertIgnoreDuplicates(opsForCvesMergedIntoMapWtRS, dbsnpCvoeClass, opsForCvesMergedIntoMapWtRSColl)
    }
}