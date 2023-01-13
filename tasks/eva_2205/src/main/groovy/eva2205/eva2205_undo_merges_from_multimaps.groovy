package eva2205

import org.springframework.data.mongodb.core.BulkOperations
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.deprecate.Application
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
def logger = LoggerFactory.getLogger(GenericApplication.class)

def ssWithMapWtMergeOps = [svoeClass, dbsnpSvoeClass].collect{new RetryableBatchingCursor<>(
        where("inactiveObjects.seq").is(options.assembly).and("eventType").is(EventType.UPDATED.toString()).and(
                "reason").regex(".* merged into rs.*").and("inactiveObjects.mapWeight").exists(true),
        dbEnv.mongoTemplate, it)}
ssWithMapWtMergeOps.each {it.each{List<SubmittedVariantOperationEntity> ops ->
    def (sveBulkUpdates, dbsnpSveBulkUpdates) = [sveClass, dbsnpSveClass].collect{
        dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, it)}
    def opsGroupedByHash = ops.groupBy{it.inactiveObjects[0].hashedMessage}
    def ssHashesThatShouldHaveMapWt = ops.collect{it.inactiveObjects[0].hashedMessage}.unique()
    def (svesThatShouldGetMapWt, dbsnpSvesThatShouldGetMapWt) = [sveClass, dbsnpSveClass].collect{
        dbEnv.mongoTemplate.find(query(where("_id").in(ssHashesThatShouldHaveMapWt).and("mapWeight").exists(false)), it)}
    def updateTargetSves = {targetSvesToUpdate, bulkUpdateObj, opsClass ->
        if (targetSvesToUpdate.size() > 0) {
            def opsToInsert = new ArrayList<>()
            targetSvesToUpdate.each {sve ->
                def correspondingOp = opsGroupedByHash.get(sve.hashedMessage)[0]
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
    }
    updateTargetSves(svesThatShouldGetMapWt, sveBulkUpdates, svoeClass)
    updateTargetSves(dbsnpSvesThatShouldGetMapWt, dbsnpSveBulkUpdates, dbsnpSvoeClass)
}}
