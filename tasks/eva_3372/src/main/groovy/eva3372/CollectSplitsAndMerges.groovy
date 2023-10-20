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
    String assembly
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv

    CollectSplitsAndMerges() {}

    CollectSplitsAndMerges(String assembly, EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv) {
        this.assembly = assembly
        this.prodEnv = prodEnv
        this.devEnv = devEnv
    }

    def run = {
        def evaAndDbsnpSveCursorsProd = [sveClass, dbsnpSveClass].collect { collectionClass ->
            new RetryableBatchingCursor<>(
                    where("seq").is(assembly).and("remappedFrom").exists(true).and("rs").exists(true),
                    prodEnv.mongoTemplate, collectionClass)
        }
        def numRecordsProcessed = 0
        evaAndDbsnpSveCursorsProd.each { cursor ->
            cursor.each { List<SubmittedVariantEntity> sves ->
                detectPendingSplits(sves)
                gatherPossibleMerges(sves)
                numRecordsProcessed += sves.size()
                logger.info("${numRecordsProcessed} SS processed so far...")
            }
        }
        logger.info("Scan through SVE and dbsnpSVE collections complete")

        detectPendingMerges()

        // Counts
        def numSplits = devEnv.mongoTemplate.count(query(where("_id").regex(/^${assembly}#/)), CveHashesWithRsId.class)
        def numMerges = devEnv.mongoTemplate.count(query(where("_id").regex(/^${assembly}#/)), RsIdsWithCveHash.class)
        logger.info("Processing for ${assembly} complete: ${numSplits} pending splits, ${numMerges} pending merges")
    }

    /**
     * Detect and store pending splits: >1 SVE with same rsID and different CVE hash
     */
    def detectPendingSplits = { List<SubmittedVariantEntity> sves ->
        // Only process rsIds that we haven't yet recorded in rsHashesWithId
        def rsIds = sves.collect { it.clusteredVariantAccession }
        def existingRecords = devEnv.mongoTemplate.find(
                query(where("_id").in(rsIds.collect { "${assembly}#${it}" })), CveHashesWithRsId.class)
        def processedRsIds = existingRecords.collect { it.rsId }
        def rsIdsToProcess = rsIds.findAll { !processedRsIds.contains(it) }

        // For these rsIds, get all SVEs and corresponding hashes
        def svesInDbWithRs = [sveClass, dbsnpSveClass].collectMany { entityClass ->
            prodEnv.mongoTemplate.find(query(where("seq").is(assembly).and("rs").in(rsIdsToProcess)), entityClass)
        }
        def svesGroupedByRs = svesInDbWithRs.groupBy { it.clusteredVariantAccession }
        def recordsToInsert = []
        svesGroupedByRs.each { rsId, svesWithRs ->
            def cveHashes = svesWithRs.collect { EVAObjectModelUtils.getClusteredVariantHash(it) }.toSet()
            if (cveHashes.size() > 1) {
                recordsToInsert.add(new CveHashesWithRsId("${assembly}#${rsId}", assembly, rsId, cveHashes))
            }
        }
        def insertedCount = 0
        if (!recordsToInsert.isEmpty()) {
            insertedCount = devEnv.bulkInsertIgnoreDuplicates(recordsToInsert, CveHashesWithRsId.class)
        }
        logger.info("Inserted ${insertedCount} split candidate documents (CveHashesWithRsId)")
    }

    /**
     * First step to detecting merges: store rsIds grouped by CVE hash
     */
    def gatherPossibleMerges = { List<SubmittedVariantEntity> sves ->
        def svesGroupedByCveHash = sves.groupBy { EVAObjectModelUtils.getClusteredVariantHash(it) }
        def cveHashes = svesGroupedByCveHash.keySet()
        def existingRecords = devEnv.mongoTemplate.find(query(where("_id").in(cveHashes.collect { "${assembly}#${it}" })), RsIdsWithCveHash.class)
        def processedCveHashes = existingRecords.collect { it.cveHash }

        // Add to existing records
        if (!existingRecords.isEmpty()) {
            def bulkOps = devEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, RsIdsWithCveHash.class)
            existingRecords.each { rsIdsWithHash ->
                def svesToProcess = svesGroupedByCveHash.get(rsIdsWithHash.cveHash)
                def rsIdsToAdd = svesToProcess.collect { it.clusteredVariantAccession }
                rsIdsWithHash.rsIds.addAll(rsIdsToAdd)
                bulkOps.updateOne(query(where("_id").is("${assembly}#${rsIdsWithHash.cveHash}")), Update.update("rsIds", rsIdsWithHash.rsIds))
            }
            bulkOps.execute()
        }

        // Create new records
        def cveHashesToProcess = cveHashes.findAll { !processedCveHashes.contains(it) }
        def recordsToInsert = []
        svesGroupedByCveHash.each { hash, svesWithHash ->
            if (cveHashesToProcess.contains(hash)) {
                def rsIds = svesWithHash.collect { it.clusteredVariantAccession }.toSet()
                recordsToInsert.add(new RsIdsWithCveHash("${assembly}#${hash}", assembly, hash, rsIds))
            }
        }
        def insertedCount = 0
        if (!recordsToInsert.isEmpty()) {
            insertedCount = devEnv.bulkInsertIgnoreDuplicates(recordsToInsert, RsIdsWithCveHash)
        }
        logger.info("Inserted ${insertedCount} merge candidate documents (RsIdsWithCveHash)")
    }

    /**
     * Second step to detecting merges: remove records for hashes without multiple rsIDs associated
     */
    def detectPendingMerges = {
        def possibleMergeCursor = new RetryableBatchingCursor<>(where("_id").regex(/^${assembly}#/), devEnv.mongoTemplate, RsIdsWithCveHash.class)
        possibleMergeCursor.each { List<RsIdsWithCveHash> batch ->
            def recordsToRemove = batch.findAll { it.rsIds.size() < 2 }
            def removed = devEnv.mongoTemplate.findAllAndRemove(query(where("_id").in(recordsToRemove.collect { "${assembly}#${it.cveHash}" })), RsIdsWithCveHash.class)
            logger.info("Removed ${removed.size()} merge candidate documents (RsIdsWithCveHash)")
        }
    }
}
