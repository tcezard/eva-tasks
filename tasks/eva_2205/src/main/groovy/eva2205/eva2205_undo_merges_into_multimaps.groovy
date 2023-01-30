package eva2205

import groovy.cli.picocli.CliBuilder
import org.springframework.data.mongodb.core.BulkOperations
import org.springframework.data.mongodb.core.query.Query
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Update

import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType

import uk.ac.ebi.eva.accession.clustering.metric.ClusteringMetric
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.batch.io.ClusteredVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.EVAObjectModelUtils
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.metrics.metric.MetricCompute

import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static eva2205.eva2205_utils.*

def cli = new CliBuilder()
cli.propertiesFile(args:1, "Properties file to use for remediation", required: true)
cli.assembly(args:1, "Assembly to remediate", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

class UndoMergesIntoMultiMaps {
    EVADatabaseEnvironment dbEnv
    MetricCompute metricCompute
    String assembly
    Long dbsnpSSIDUpperBound, dbsnpRSIDUpperBound
    static def logger = LoggerFactory.getLogger(GenericApplication.class)

    UndoMergesIntoMultiMaps(EVADatabaseEnvironment dbEnv, String assembly) {
        this.dbEnv = dbEnv
        this.assembly = assembly
        this.dbsnpSSIDUpperBound = dbEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class)
        this.dbsnpRSIDUpperBound = dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class)
        this.metricCompute = dbEnv.springApplicationContext.getBean(MetricCompute.class)
    }

    /**
     * Undo merges into multimap RS
     * See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=1814908308
     * @param recentSvoeOpsWithoutMapWt - Most recent SVOE/dbsnpSVOE records without map weight for a set of SS hashes (keyed by SS hash)
     * @param currentSSRecords - Current SS record for a set of SS hashes (keyed by SS hash)
     * @return None
     */
    void undoMultiMapRSAssignment(Map<String, List> recentSvoeOpsWithoutMapWt, Map<String, List> currentSSRecords) {
        def (bulkSVEUpdateOps, bulkdbSnpSVEUpdateOps) =
        [sveClass, dbsnpSveClass].collect { dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, it) }
        def (svoeInserts, dbSnpSvoeInserts, dbsnpCvoeInserts) = [new ArrayList<>(), new ArrayList<>(),
                                                                 new ArrayList<>()]
        def (numSVEUpdates, numDbsnpSveUpdates) = [0, 0]

        recentSvoeOpsWithoutMapWt.each { ssHash, recentSSOpWithoutMapWt ->
            def ssRecordToUpdate = currentSSRecords.get(ssHash)[0]
            def oldRS = ssRecordToUpdate.clusteredVariantAccession
            // Since we are undoing a previous RS assignment to this SS,
            // the new RS is the one from the Op recording that previous assignment
            def newRS = recentSSOpWithoutMapWt.inactiveObjects[0].clusteredVariantAccession
            Query queryToFindSSToUpdate = query(where("_id").is(ssHash).and("rs").ne(newRS))

            def ssAccession = ssRecordToUpdate.accession
            def isEvaSS = (ssAccession >= dbsnpSSIDUpperBound)
            def (sveUpdateObj, svoeInsertObj) = isEvaSS ?
                    [bulkSVEUpdateOps, svoeInserts] : [bulkdbSnpSVEUpdateOps, dbSnpSvoeInserts]

            logger.info("SS with hash ${ssHash} and accession ${ssAccession} will be assigned " +
                    "non-map weighted RS ${newRS}...")
            def updatesToMake = new Update()
            updatesToMake.set("rs", newRS)
            updatesToMake.unset("mapWeight")
            sveUpdateObj.updateOne(queryToFindSSToUpdate, updatesToMake)

            def svoeOpToRecordUnmerge = new SubmittedVariantOperationEntity()
            svoeOpToRecordUnmerge.fill(EventType.UPDATED, ssAccession, "Undo map-weighted RS assignment.",
                    Arrays.asList(new SubmittedVariantInactiveEntity(ssRecordToUpdate)))
            svoeOpToRecordUnmerge.setId("UNDO_MULTIMAP_RS_ASSIGN_${ssRecordToUpdate.hashedMessage}")
            svoeInsertObj.add(svoeOpToRecordUnmerge)
            if (isEvaSS) {
                numSVEUpdates += 1
            } else {
                numDbsnpSveUpdates += 1
            }

            def cvoeOpToRecordUnmerge = new ClusteredVariantOperationEntity()
            def oldCVE = EVAObjectModelUtils.toClusteredVariantEntity(ssRecordToUpdate)
            cvoeOpToRecordUnmerge.fill(EventType.UNDO_MERGE, oldRS, newRS, "Undo merge into map-weighted RS.",
                    Arrays.asList(new ClusteredVariantInactiveEntity(oldCVE)))
            cvoeOpToRecordUnmerge.setId("UNDO_MERGE_INTO_MULTIMAP_${oldCVE.hashedMessage}")
            dbsnpCvoeInserts.add(cvoeOpToRecordUnmerge)
        }
        if (numSVEUpdates > 0) bulkSVEUpdateOps.execute()
        if (numDbsnpSveUpdates > 0) bulkdbSnpSVEUpdateOps.execute()

        dbEnv.bulkInsertIgnoreDuplicates(svoeInserts, svoeClass)
        dbEnv.bulkInsertIgnoreDuplicates(dbSnpSvoeInserts, dbsnpSvoeClass)
        dbEnv.bulkInsertIgnoreDuplicates(dbsnpCvoeInserts, dbsnpCvoeClass)
    }

    List undoSVEsAssignedMultimaps(List mapWtSVEs) {
        def ssHashesWithMapWt = mapWtSVEs.collect { it.hashedMessage }
        // See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=2145326284
        def ssOpsInvolvedInMergeIntoMultimaps = [svoeClass, dbsnpSvoeClass].collect {
            dbEnv.mongoTemplate.find(query(where("inactiveObjects.seq").is(assembly)
                    .and("eventType").is(EventType.UPDATED)
                    .and("reason").regex(".* was merged into .*")
                    .and("inactiveObjects.hashedMessage").in(ssHashesWithMapWt)), it)}.flatten()
        def impactedSSOpHistory = ssOpsInvolvedInMergeIntoMultimaps.groupBy { it.inactiveObjects[0].hashedMessage }
        def recentSSRecordsWithoutMapWt = impactedSSOpHistory.collectEntries { ssHash, opsForSSHash ->
            def mergeChain = getChronologicalMergeChain(opsForSSHash)
            if (mergeChain.size() > 1) {
                logger.info("These operations involve more than one merge in a merge chain:" + mergeChain)
            }
            def originalSSRecord = mergeChain[0]
            // Return the most recent operation that recorded the SS without a map-weight,
            // if the original SS did not have a map-weight
            if (Objects.nonNull(originalSSRecord) && Objects.isNull(originalSSRecord.inactiveObjects[0].mapWeight)) {
                return [ssHash, mostRecentNonMapWtSSRecord(mergeChain)]
            }
            return [ssHash, null]
        }.findAll { k, v -> Objects.nonNull(v) }
        undoMultiMapRSAssignment(recentSSRecordsWithoutMapWt, mapWtSVEs.groupBy{it.hashedMessage})
        return recentSSRecordsWithoutMapWt.values().flatten()
    }

    def deprecateAndClearSSBatch(ssToDeprecate) {
        DeprecateMapWtSS.deprecateSS(dbEnv, ssToDeprecate.toSet()) // Use unique to weed out duplicates within a batch
        ssToDeprecate.clear()
    }

    def deprecateOrphanedCvesWithSameHash(cvesToResurrect) {
        def cvesWithSameHash = [cveClass, dbsnpCveClass].collect{coll ->
            dbEnv.mongoTemplate.find(query(where("_id").in(cvesToResurrect.collect{it.hashedMessage}).and(
                    "mapWeight").exists(true)), coll)
        }.flatten()
        def cvDeprecationWriter = new ClusteredVariantDeprecationWriter(assembly,
                dbEnv.mongoTemplate,
                dbEnv.submittedVariantAccessioningService,
                dbsnpRSIDUpperBound,
                "EVA2205", "Variant deprecated due to mapWeight > 1")
        // By the design of the write, this will only succeed if these RS IDs are orphaned i.e., there are no SS that have these RS IDs
        cvDeprecationWriter.write(cvesWithSameHash)
        this.metricCompute.addCount(ClusteringMetric.CLUSTERED_VARIANTS_DEPRECATED,
                cvDeprecationWriter.numDeprecatedEntities)
        this.metricCompute.saveMetricsCountsInDB()
    }

    // See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=1814908308
    def undoMergesIntoMultiMapRS() {
        def mapWtSSCursors = [sveClass, dbsnpSveClass].collect {
            new RetryableBatchingCursor(where("seq").is(assembly).and("mapWeight").exists(true), dbEnv.mongoTemplate, it)
        }
        def collectionWithCvesToResurrect = "cvesToResurrect"
        def svesStillWithMapWtRS = new ArrayList<>()
        def batchIndex = 0
        mapWtSSCursors.each {it.each {mapWtSVEs ->
            def involvedRS = mapWtSVEs.findAll{
                Objects.nonNull(it.clusteredVariantAccession)}.collect{it.clusteredVariantAccession}.toSet()
            def involvedHashes = mapWtSVEs.collect{it.hashedMessage}
            // Also pull in other SS with the same RS so that both SS and RS deprecation can take place within a batch
            mapWtSVEs += [sveClass, dbsnpSveClass].collect{collectionClass ->
                dbEnv.mongoTemplate.find(query(where("_id").nin(involvedHashes).and("seq").is(assembly)
                        .and("mapWeight").exists(true).and("rs").in(involvedRS)), collectionClass)}.flatten()
            def svesWithUndoneMultimaps = undoSVEsAssignedMultimaps(mapWtSVEs)
            def sveHashesWithUndoneMultimaps = svesWithUndoneMultimaps.collect{
                it.inactiveObjects[0].hashedMessage}.toSet()
            // Restore CVEs which were incorrectly merged into map weighted CVEs
            def cvesToResurrect = svesWithUndoneMultimaps.collect {
                EVAObjectModelUtils.toClusteredVariantEntity(it.inactiveObjects[0].toSubmittedVariantEntity())
            }
            // Store these separately because these cannot be obtained again with a subsequent run in case of failures.
            // This is because the deprecation below will take out the SVEs with mapWeight entries and that will impact mapWtSSCursor above.
            dbEnv.bulkInsertIgnoreDuplicates(cvesToResurrect, cveClass, collectionWithCvesToResurrect)
            // All the legitimate non map-weight RS assignments to SS have been carried out by getCVEsToResurrect above
            // Therefore, we can proceed to make room for resurrected RS by deprecating SS associated with the current set of map-weighted RS
            svesStillWithMapWtRS = mapWtSVEs.findAll{!sveHashesWithUndoneMultimaps.contains(it.hashedMessage)}
            logger.info("Deprecating ${svesStillWithMapWtRS.size()} SS in batch ${batchIndex}...")
            deprecateAndClearSSBatch(svesStillWithMapWtRS)
            batchIndex += 1
        }}

        // Resurrect CVEs that were collected above
        def cvesToResurrectCursor = new RetryableBatchingCursor(where("asm").is(assembly), dbEnv.mongoTemplate,
                dbsnpCveClass, 1000, collectionWithCvesToResurrect)
        cvesToResurrectCursor.each {cvesToResurrect ->
            // It is possible that all constituent SS of a map-weighted RS were assigned non map-wt RS above in undoSVEsAssignedMultimaps
            // In that case, the SS deprecation above won't take place (because only remaining map-weighted SS are chosen for deprecation)
            // and hence the older map-weighted RS assigned to these SS will be orphaned and won't be deprecated either.
            // Explicitly deprecate such RS.
            deprecateOrphanedCvesWithSameHash(cvesToResurrect)
            def numCvesResurrected = dbEnv.bulkInsertIgnoreDuplicates(cvesToResurrect.findAll {
                it.accession < dbsnpRSIDUpperBound }, dbsnpCveClass)
            numCvesResurrected += dbEnv.bulkInsertIgnoreDuplicates(cvesToResurrect.findAll
                    { it.accession >= dbsnpRSIDUpperBound }, cveClass)
            this.metricCompute.addCount(ClusteringMetric.CLUSTERED_VARIANTS_CREATED, numCvesResurrected)
        }
        logger.info("Resurrected ${this.metricCompute.getCount(ClusteringMetric.CLUSTERED_VARIANTS_CREATED)} CVEs" +
                "for assembly ${this.assembly}...")
        this.metricCompute.saveMetricsCountsInDB()
    }
}

EVADatabaseEnvironment dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assembly])
def undoMergeObj = new UndoMergesIntoMultiMaps(dbEnv, options.assembly)
undoMergeObj.undoMergesIntoMultiMapRS()
