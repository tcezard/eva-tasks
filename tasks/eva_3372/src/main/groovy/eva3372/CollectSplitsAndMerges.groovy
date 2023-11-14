package eva3372

import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.BulkOperations
import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.getDbsnpSveClass
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.getSveClass

class CollectSplitsAndMerges {
    static def logger = LoggerFactory.getLogger(CollectSplitsAndMerges.class)

    public static final String PENDING_MERGES_COLLECTION_NAME = "pendingMerges"
    public static final String PENDING_SPLITS_COLLECTION_NAME = "pendingSplits"

    String assembly
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv

    CollectSplitsAndMerges() {}

    CollectSplitsAndMerges(String assembly, EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv) {
        this.assembly = assembly
        this.prodEnv = prodEnv
        this.devEnv = devEnv
    }

    def collectAndCount = {
        def evaAndDbsnpSveCursorsProd = [sveClass, dbsnpSveClass].collect { collectionClass ->
            new RetryableBatchingCursor<>(
                    where("seq").is(assembly).and("rs").exists(true),
                    prodEnv.mongoTemplate, collectionClass)
        }
        def numRecordsProcessed = 0
        evaAndDbsnpSveCursorsProd.each { cursor ->
            cursor.each { List<SubmittedVariantEntity> sves ->
                def remappedSves = sves.findAll { !it.remappedFrom.equals("") }
                if (!remappedSves.isEmpty()) {
                    gatherPossibleSplits(remappedSves)
                    gatherPossibleMerges(remappedSves)
                    numRecordsProcessed += remappedSves.size()
                    logger.info("${numRecordsProcessed} SS processed so far...")
                }
            }
        }
        logger.info("Scan through SVE and dbsnpSVE collections complete")

        detectPendingSplits()
        detectPendingMerges()

        // Counts
        def numSplits = devEnv.mongoTemplate.count(query(where("_id").regex(/^${assembly}#/)), CveHashesWithRsId.class, PENDING_SPLITS_COLLECTION_NAME)
        def numMerges = devEnv.mongoTemplate.count(query(where("_id").regex(/^${assembly}#/)), RsIdsWithCveHash.class, PENDING_MERGES_COLLECTION_NAME)
        logger.info("Processing for ${assembly} complete: ${numSplits} pending splits, ${numMerges} pending merges")
    }

    /**
     * First step to detecting splits: store CVE hashes grouped by RS ID
     */
    def gatherPossibleSplits = { List<SubmittedVariantEntity> sves ->
        def svesGroupedByRsId = sves.groupBy { it.clusteredVariantAccession }
        def rsIds = svesGroupedByRsId.keySet()
        def existingRecords = devEnv.mongoTemplate.find(
                query(where("_id").in(rsIds.collect { "${assembly}#${it}".toString() })), CveHashesWithRsId.class)

        // Add to existing records
        if (!existingRecords.isEmpty()) {
            boolean needsUpdate = false
            def bulkOps = devEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, CveHashesWithRsId.class)
            existingRecords.each { cveHashesWithRsId ->
                def svesToProcess = svesGroupedByRsId.get(cveHashesWithRsId.rsId)
                def cveHashesToAdd = svesToProcess.collect { EVAObjectModelUtils.getClusteredVariantHash(it) }
                if (!cveHashesWithRsId.cveHashes.containsAll(cveHashesToAdd)) {
                    cveHashesWithRsId.cveHashes.addAll(cveHashesToAdd)
                    bulkOps.updateOne(query(where("_id").is(cveHashesWithRsId.id)), Update.update("cveHashes", cveHashesWithRsId.cveHashes))
                    needsUpdate = true
                }
            }
            if (needsUpdate) {
                def bulkResult = bulkOps.execute()
                if (bulkResult.modifiedCountAvailable && bulkResult.modifiedCount > 0) {
                    logger.info("Updated ${bulkResult.modifiedCount} split candidate documents (CveHashesWithRsId)")
                }
            }
        }

        // Create new records
        def processedRsIds = existingRecords.collect { it.rsId }
        def rsIdsToProcess = rsIds.findAll { !processedRsIds.contains(it) }
        def recordsToInsert = []
        svesGroupedByRsId.each { rsId, svesWithRsId ->
            if (rsIdsToProcess.contains(rsId)) {
                def cveHashes = svesWithRsId.collect { EVAObjectModelUtils.getClusteredVariantHash(it) }.toSet()
                recordsToInsert.add(new CveHashesWithRsId("${assembly}#${rsId}", assembly, rsId, cveHashes))
            }
        }
        if (!recordsToInsert.isEmpty()) {
            def insertedCount = devEnv.bulkInsertIgnoreDuplicates(recordsToInsert, CveHashesWithRsId.class)
            logger.info("Inserted ${insertedCount} split candidate documents (CveHashesWithRsId)")
        }

    }

    /**
     * First step to detecting merges: store rsIds grouped by CVE hash
     */
    def gatherPossibleMerges = { List<SubmittedVariantEntity> sves ->
        def svesGroupedByCveHash = sves.groupBy { EVAObjectModelUtils.getClusteredVariantHash(it) }
        def cveHashes = svesGroupedByCveHash.keySet()
        def existingRecords = devEnv.mongoTemplate.find(query(where("_id").in(cveHashes.collect { "${assembly}#${it}".toString() })), RsIdsWithCveHash.class)

        // Add to existing records
        if (!existingRecords.isEmpty()) {
            boolean needsUpdate = false
            def bulkOps = devEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, RsIdsWithCveHash.class)
            existingRecords.each { rsIdsWithHash ->
                def svesToProcess = svesGroupedByCveHash.get(rsIdsWithHash.cveHash)
                def rsIdsToAdd = svesToProcess.collect { it.clusteredVariantAccession }
                if (!rsIdsWithHash.rsIds.containsAll(rsIdsToAdd)) {
                    rsIdsWithHash.rsIds.addAll(rsIdsToAdd)
                    bulkOps.updateOne(query(where("_id").is(rsIdsWithHash.id)), Update.update("rsIds", rsIdsWithHash.rsIds))
                    needsUpdate = true
                }
            }
            if (needsUpdate) {
                def bulkResult = bulkOps.execute()
                if (bulkResult.modifiedCountAvailable && bulkResult.modifiedCount > 0) {
                    logger.info("Updated ${bulkResult.modifiedCount} merge candidate documents (RsIdsWithCveHash)")
                }
            }
        }

        // Create new records
        def processedCveHashes = existingRecords.collect { it.cveHash }
        def cveHashesToProcess = cveHashes.findAll { !processedCveHashes.contains(it) }
        def recordsToInsert = []
        svesGroupedByCveHash.each { hash, svesWithHash ->
            if (cveHashesToProcess.contains(hash)) {
                def rsIds = svesWithHash.collect { it.clusteredVariantAccession }.toSet()
                recordsToInsert.add(new RsIdsWithCveHash("${assembly}#${hash}", assembly, hash, rsIds))
            }
        }
        if (!recordsToInsert.isEmpty()) {
            def insertedCount = devEnv.bulkInsertIgnoreDuplicates(recordsToInsert, RsIdsWithCveHash.class)
            logger.info("Inserted ${insertedCount} merge candidate documents (RsIdsWithCveHash)")
        }
    }

    /**
     * Second step to detecting splits: remove records for rsIDs without multiple CVE hashes associated
     */
    def detectPendingSplits = {
        def possibleSplitCursor = new RetryableBatchingCursor<>(where("_id").regex(/^${assembly}#/), devEnv.mongoTemplate, CveHashesWithRsId.class)
        possibleSplitCursor.each { List<CveHashesWithRsId> batch ->
            def recordsToCopy = batch.findAll { it.cveHashes.size() > 1 }
            if (!recordsToCopy.isEmpty()) {
                def copiedCount = devEnv.bulkInsertIgnoreDuplicates(recordsToCopy, CveHashesWithRsId.class, PENDING_SPLITS_COLLECTION_NAME)
                logger.info("Copied ${copiedCount} split candidate documents (CveHashesWithRsId)")
            }
        }
        devEnv.mongoTemplate.dropCollection("cveHashesWithRsId")
    }

    /**
     * Second step to detecting merges: remove records for hashes without multiple rsIDs associated
     */
    def detectPendingMerges = {
        def possibleMergeCursor = new RetryableBatchingCursor<>(where("_id").regex(/^${assembly}#/), devEnv.mongoTemplate, RsIdsWithCveHash.class)
        possibleMergeCursor.each { List<RsIdsWithCveHash> batch ->
            def recordsToCopy = batch.findAll { it.rsIds.size() > 1 }
            if (!recordsToCopy.isEmpty()) {
                def copiedCount = devEnv.bulkInsertIgnoreDuplicates(recordsToCopy, RsIdsWithCveHash.class, PENDING_MERGES_COLLECTION_NAME)
                logger.info("Copied ${copiedCount} merge candidate documents (RsIdsWithCveHash)")
            }
        }
        devEnv.mongoTemplate.dropCollection("rsIdsWithCveHash")
    }
}
