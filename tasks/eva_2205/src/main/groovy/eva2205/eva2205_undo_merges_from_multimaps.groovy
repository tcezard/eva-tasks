package eva2205

import org.springframework.data.mongodb.core.BulkOperations
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static eva2205.eva2205_utils.*
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


class UndoMergesFromMultiMaps {
    EVADatabaseEnvironment dbEnv
    String assembly
    static def logger = LoggerFactory.getLogger(GenericApplication.class)

    UndoMergesFromMultiMaps(EVADatabaseEnvironment dbEnv, String assembly) {
        this.dbEnv = dbEnv
        this.assembly = assembly
    }

    private void updateTargetSves(ssHashesAndMapWtOps, targetSvesToUpdate, bulkUpdateObj, opsClass) {
        def opsToInsert = new ArrayList<>()
        targetSvesToUpdate.each { sve ->
            def correspondingOp = ssHashesAndMapWtOps.get(sve.hashedMessage)
            def mapWeightedRSToAssign = correspondingOp.inactiveObjects[0].clusteredVariantAccession
            def updatesToMake = new Update()
            updatesToMake.set("mapWeight", correspondingOp.inactiveObjects[0].mapWeight)
            updatesToMake.set("rs", mapWeightedRSToAssign)
            bulkUpdateObj.updateOne(query(where("_id").is(sve.hashedMessage)), updatesToMake)
            logger.info("SS with hash ${sve.hashedMessage} and accession ${sve.accession} will be assigned " +
                    "map-weighted RS ${mapWeightedRSToAssign}")

            def opToRecordRestoredMapWt = new SubmittedVariantOperationEntity()
            opToRecordRestoredMapWt.fill(EventType.UPDATED, sve.accession, "Restore map-weighted RS ID.",
                    Arrays.asList(new SubmittedVariantInactiveEntity(sve)))
            opToRecordRestoredMapWt.setId("RESTORE_MAPWT_${sve.hashedMessage}")
            opsToInsert.add(opToRecordRestoredMapWt)
        }
        bulkUpdateObj.execute()
        dbEnv.bulkInsertIgnoreDuplicates(opsToInsert, opsClass)
    }

    void undoMergesFromMultimapRS() {
        def ssWithMapWtMergeOps = [svoeClass, dbsnpSvoeClass].collect {
            new RetryableBatchingCursor<>(
                    where("inactiveObjects.seq").is(assembly).and("eventType").is(EventType.UPDATED.toString()).and(
                            "reason").regex(".* merged into rs.*").and("inactiveObjects.mapWeight").exists(true),
                    dbEnv.mongoTemplate, it)
        }
        ssWithMapWtMergeOps.each {it.each { ops ->
            def ssHashesAndMapWtOps = getSSHistoryInvolvedInRSMerges(dbEnv, assembly,ops)
                    .groupBy { it.inactiveObjects[0].hashedMessage }.collectEntries { ssHash, opsForSSHash ->
                def mergeChain = getChronologicalMergeChain(opsForSSHash)
                if (mergeChain.size() > 1) {
                    logger.info("These operations involve more than one merge in a merge chain:" + mergeChain)
                }
                def originalSSRecord = mergeChain[0]
                // Return the most recent operation that recorded the SS with a map-weight, if the original SS has a map-weight
                if (Objects.nonNull(originalSSRecord) && Objects.nonNull(originalSSRecord.inactiveObjects[0].mapWeight)) {
                    return [ssHash, mostRecentMapWtSSRecord(mergeChain)]
                }
                return [ssHash, null]
            }.findAll { k, v -> Objects.nonNull(v) }

            def (sveBulkUpdates, dbsnpSveBulkUpdates) = [sveClass, dbsnpSveClass].collect {
                dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, it)
            }
            def ssHashesThatShouldHaveMapWt = ssHashesAndMapWtOps.keySet()
            def (svesThatShouldGetMapWt, dbsnpSvesThatShouldGetMapWt) = [sveClass, dbsnpSveClass].collect {
                dbEnv.mongoTemplate.find(query(where("_id").in(ssHashesThatShouldHaveMapWt)
                        .and("mapWeight").exists(false)), it)
            }

            if (svesThatShouldGetMapWt.size() > 0) {
                updateTargetSves(ssHashesAndMapWtOps, svesThatShouldGetMapWt, sveBulkUpdates, svoeClass)
            }
            if (dbsnpSvesThatShouldGetMapWt.size() > 0) {
                updateTargetSves(ssHashesAndMapWtOps, dbsnpSvesThatShouldGetMapWt, dbsnpSveBulkUpdates, dbsnpSvoeClass)
            }
        }}
    }
}

EVADatabaseEnvironment dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assembly])
def undoMergeObj = new UndoMergesFromMultiMaps(dbEnv, options.assembly)
undoMergeObj.undoMergesFromMultimapRS()
