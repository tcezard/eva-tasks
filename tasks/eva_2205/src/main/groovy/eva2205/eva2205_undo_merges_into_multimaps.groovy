package eva2205

import eva2205.DeprecateMapWtSS
import org.springframework.data.mongodb.core.BulkOperations
import org.springframework.data.mongodb.core.query.Criteria
import org.springframework.data.mongodb.core.query.Query
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVAObjectModelUtils
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import groovy.cli.picocli.CliBuilder

def cli = new CliBuilder()
cli.propertiesFile(args:1, "Properties file to use for remediation", required: true)
cli.assembly(args:1, "Assembly to remediate", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assembly])
def dbsnpSSIDUpperBound = dbEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class)
def dbsnpRSIDUpperBound = dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class)
def logger = LoggerFactory.getLogger(GenericApplication.class)

def undoMergeIntoMultiMapRS = {svoeOps ->
    def (bulkSVEUpdateOps, bulkdbSnpSVEUpdateOps) =
    [sveClass, dbsnpSveClass].collect { dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, it) }
    def (svoeInserts, dbSnpSvoeInserts) = [new ArrayList<>(), new ArrayList<>()]
    def (numSVEUpdates, numDbsnpSveUpdates, numDbsnpCvoeUpdates) = [0, 0, 0]

    svoeOps.each {SubmittedVariantOperationEntity svoeOp ->
        def ssHash = svoeOp.inactiveObjects[0].hashedMessage
        Query queryToFindSSToUpdate = query(where("_id").is(ssHash).and("rs").exists(true))
        def ssRecordsToUpdate = [sveClass, dbsnpSveClass].collect{
            dbEnv.mongoTemplate.find(queryToFindSSToUpdate, it)}.flatten()
        if (ssRecordsToUpdate.size() > 0) {
            def ssAccession = ssRecordsToUpdate[0].accession
            def isEvaSS = (ssAccession >= dbsnpSSIDUpperBound)
            def (sveUpdateObj, svoeInsertObj) =  isEvaSS?
                    [bulkSVEUpdateOps, svoeInserts] :[bulkdbSnpSVEUpdateOps, dbSnpSvoeInserts]
            def oldRS = ssRecordsToUpdate[0].clusteredVariantAccession
            // Since we are undoing a previous RS assignment to this SS,
            // the new RS is the one from the Op recording that previous assignment
            def newRS = svoeOp.inactiveObjects[0].clusteredVariantAccession

            logger.info("SS with hash ${ssHash} and accession ${ssAccession} will be assigned " +
                    "non-map weighted RS ${newRS}...")
            sveUpdateObj.updateOne(queryToFindSSToUpdate, Update.update("rs", newRS))
            sveUpdateObj.updateOne(queryToFindSSToUpdate, new Update().unset("mapWeight"))

            def svoeOpToRecordUnmerge = new SubmittedVariantOperationEntity()
            svoeOpToRecordUnmerge.fill(EventType.UPDATED, ssAccession, "Undo map-weighted RS assignment.",
                    Arrays.asList(new SubmittedVariantInactiveEntity(ssRecordsToUpdate[0])))
            svoeOpToRecordUnmerge.setId("UNDO_MULTIMAP_RS_ASSIGN_${ssRecordsToUpdate[0].hashedMessage}")
            svoeInsertObj.add(svoeOpToRecordUnmerge)
            if(isEvaSS) { numSVEUpdates += 1 } else { numDbsnpSveUpdates += 1 }

            def cvoeOpToRecordUnmerge = new ClusteredVariantOperationEntity()
            def oldCVE= EVAObjectModelUtils.toClusteredVariantEntity(ssRecordsToUpdate[0])
            cvoeOpToRecordUnmerge.fill(EventType.MERGED, oldRS, newRS, "Undo merge into map-weighted RS.",
                    Arrays.asList(new ClusteredVariantInactiveEntity(oldCVE)))
            cvoeOpToRecordUnmerge.setId("UNDO_MERGE_INTO_MULTIMAP_${oldCVE.hashedMessage}")
            dbsnpCvoeInserts.add(cvoeOpToRecordUnmerge)
        }
    }
    if(numSVEUpdates > 0) bulkSVEUpdateOps.execute()
    if(numDbsnpSveUpdates > 0) bulkdbSnpSVEUpdateOps.execute()

    dbEnv.bulkInsertIgnoreDuplicates(svoeInserts, svoeClass)
    dbEnv.bulkInsertIgnoreDuplicates(dbSnpSvoeInserts, dbsnpSvoeClass)
    dbEnv.bulkInsertIgnoreDuplicates(dbsnpCvoeInserts, dbsnpCvoeClass)
}

def getCVEsToResurrect = {mapWtCVEs ->

    def mapWtRSIDs = mapWtCVEs.collect{it.accession}
    // See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=2145326284
    def allCVEsMergedIntoMapWtCVEs = [cvoeClass, dbsnpCvoeClass].collect{dbEnv.mongoTemplate.find(query(
            where("inactiveObjects.asm").is(options.assembly).and("eventType").is(EventType.MERGED).and(
                    "mergeInto").in(mapWtRSIDs)), it)}.flatten()
    def nonMapWtCVEsMergedIntoMapWtCVEs = allCVEsMergedIntoMapWtCVEs.findAll{ Objects.isNull(it.inactiveObjects[0].mapWeight)}
    def nonMapWtRSIDs = nonMapWtCVEsMergedIntoMapWtCVEs.collect{it.accession}
    def reasonRegEx = mapWtRSIDs.collect{".* merged into rs${it}."}.join("|")
    def nonMapWtSSAssignedMapWtRS = [svoeClass, dbsnpSvoeClass].collect{dbEnv.mongoTemplate.find(query(where("inactiveObjects.seq").is(
            options.assembly).and("inactiveObjects.rs").in(nonMapWtRSIDs).and("eventType").is(EventType.UPDATED).and(
            "reason").regex(reasonRegEx).and("inactiveObjects.mapWeight").exists(false)), it)}.flatten()
    // Restore CVEs which were incorrectly merged into map weighted CVEs
    def cvesToResurrect = nonMapWtCVEsMergedIntoMapWtCVEs.collect{it.inactiveObjects[0].toClusteredVariantEntity()}
    undoMergeIntoMultiMapRS(nonMapWtSSAssignedMapWtRS)

    // If there were other multimap RS merged into multimap RS in the current batch
    // recursively follow these ancestral merges to find any non-multimap RS in the past that should be resurrected
    def mapWtCVEsMergedIntoMapWtRS = allCVEsMergedIntoMapWtCVEs.findAll{ Objects.nonNull(
            it.inactiveObjects[0].mapWeight)}.collect{it.inactiveObjects[0].toClusteredVariantEntity()}
    if (mapWtCVEsMergedIntoMapWtRS.size() > 0) {
        cvesToResurrect.addAll(getCVEsToResurrect(mapWtCVEsMergedIntoMapWtRS))
    }
    return cvesToResurrect
}

def mapWtRSCursor = new RetryableBatchingCursor(where("asm").is(options.assembly).and("mapWeight").exists(true),
        dbEnv.mongoTemplate, dbsnpCveClass)
def collectionWithCvesToResurrect = "cvesToResurrect"
mapWtRSCursor.each {mapWtCVEs ->
    def cvesToResurrect = getCVEsToResurrect(mapWtCVEs)
    // Store these separately because these cannot be obtained again with a subsequent run in case of failures.
    // This is because the deprecation below will take out the CVEs with mapWeight entries and that will impact mapWtRSCursor above.
    dbEnv.bulkInsertIgnoreDuplicates(cvesToResurrect, cveClass, collectionWithCvesToResurrect)
    // All the legitimate non map-weight RS assignments to SS have been carried out by getCVEsToResurrect above
    // Therefore, we can proceed to make room for resurrected RS by deprecating SS associated with the current set of map-weighted RS
    def svesWithMapWtRS = [sveClass, dbsnpSveClass].collect{dbEnv.mongoTemplate.find(query(where("seq").is(options.assembly).and("mapWeight").exists(true).and(
            "rs").in(mapWtCVEs.collect{it.accession})), it)}.flatten()
    DeprecateMapWtSS.deprecateSS(dbEnv, svesWithMapWtRS)
}

// Resurrect CVEs that were collected above
def cvesToResurrectCursor = new RetryableBatchingCursor(new Criteria(), dbEnv.mongoTemplate, dbsnpCveClass, 1000,
        collectionWithCvesToResurrect)
cvesToResurrectCursor.each {cvesToResurrect ->
    dbEnv.bulkInsertIgnoreDuplicates(cvesToResurrect.findAll { it.accession < dbsnpRSIDUpperBound }, dbsnpCveClass)
    dbEnv.bulkInsertIgnoreDuplicates(cvesToResurrect.findAll { it.accession >= dbsnpRSIDUpperBound }, cveClass)
}
