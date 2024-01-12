package eva3450

import org.slf4j.LoggerFactory
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.clustering.configuration.batch.io.RSMergeAndSplitCandidatesReaderConfiguration
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.getDbsnpSveClass
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.getSveClass

class CreateSplitCandidates {
    static def logger = LoggerFactory.getLogger(CreateSplitCandidates.class)

    String assembly
    String splitCandidatesFile
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv
    int batchSize

    CreateSplitCandidates() {}

    CreateSplitCandidates(String assembly, String splitCandidatesFile, EVADatabaseEnvironment prodEnv,
                          EVADatabaseEnvironment devEnv, int batchSize) {
        this.assembly = assembly
        this.splitCandidatesFile = splitCandidatesFile
        this.prodEnv = prodEnv
        this.devEnv = devEnv
        this.batchSize = batchSize
    }


    def process = {
        def is = new File(splitCandidatesFile).newInputStream()
        logger.info("Create Split candidate documents in batches of ${batchSize}")
        int totalInsertedCount = 0
        int lineRead = 0
        def splitCandidateOperations = []
        ArrayList<String> batchedLines = new ArrayList<String>()
        is.eachLine{
            lineRead = lineRead + 1
            batchedLines.add(it)
            if (batchedLines.size() == batchSize){
                splitCandidateOperations = createBatchSplitCandidates(batchedLines)
                totalInsertedCount = insertOperations(splitCandidateOperations, totalInsertedCount, lineRead)
                batchedLines = new ArrayList()
            }
        }
        if (batchedLines.size() > 0){
            splitCandidateOperations = createBatchSplitCandidates(batchedLines)
            totalInsertedCount = insertOperations(splitCandidateOperations, totalInsertedCount, lineRead)
        }
        is.close()

    }

    def insertOperations = {ArrayList splitCandidateOperations, int totalInsertedCount, int lineRead ->
        def insertedCount = devEnv.bulkInsertIgnoreDuplicates(splitCandidateOperations, SubmittedVariantOperationEntity.class)
        totalInsertedCount = totalInsertedCount + insertedCount
        logger.info("Read ${lineRead}, Inserted ${totalInsertedCount} Split candidate documents")
        return totalInsertedCount
    }

    def createBatchSplitCandidates = { List<String> rsidAndHashes ->
        Map<Long, String[]> rsidsToHashes = rsidAndHashes.collectEntries {
            String[] sp_line = it.split("\\t")
            def rsid = sp_line.first().toLong()
            def hashes = sp_line.tail()
            return [ (rsid): hashes]
        }
        def rsids = rsidsToHashes.keySet().toList()
        def evaAndDbsnpSveCursors = [sveClass, dbsnpSveClass].collect { collectionClass ->
            prodEnv.mongoTemplate.find(query(where("seq").is(assembly).and("rs").in(rsids)), collectionClass)
        }
        def allSves = evaAndDbsnpSveCursors.each { cursor ->
            cursor.each { sves -> sves}}.flatten()
        Map<Long, SubmittedVariantEntity[]> rsidToSves = new HashMap<>()
        allSves.each {SubmittedVariantEntity sve ->
            if (! rsidToSves.containsKey(sve.clusteredVariantAccession) ){
                rsidToSves.put(sve.clusteredVariantAccession, new ArrayList<SubmittedVariantEntity>())
            }
            rsidToSves[sve.clusteredVariantAccession].add(sve)
        }
        def splitCandidateOperations= rsids.collect {Long rsid ->
            def sves = rsidToSves.get(rsid)
            def hashes = rsidsToHashes.get(rsid)
            if (! isSplitIsValid(sves, hashes, rsid)) {
                return null
            }
            // create SPLIT_CANDIDATE operation
            def submittedVariantInactiveEntity = allSves.collect {new SubmittedVariantInactiveEntity(it)}
            SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity()
            svoe.fill(
                    EventType.RS_SPLIT_CANDIDATES, rsid,
                    "Hash mismatch with " + rsid, submittedVariantInactiveEntity
            )
            svoe.setId(String.format("%s_%s_%s", RSMergeAndSplitCandidatesReaderConfiguration.SPLIT_CANDIDATE_ID_PREFIX,
                    assembly, rsid))
            return svoe
        }
        // Remove null for invalid SPLIT candidates
        splitCandidateOperations.removeAll([null])
        return splitCandidateOperations
    }

    def isSplitIsValid = { ArrayList<SubmittedVariantEntity> allSves, String[] hashes, Long rsid ->
        //Check that all the sve provided have the same clustered variant hash which match the one from the input
        def hashSet = allSves.collect { sve -> EVAObjectModelUtils.getClusteredVariantHash(sve) }.toSet()
        if (hashSet == hashes.toList().toSet()){
            return true
        } else {
            logger.warn("Hash from files ${hashSet} is different from Hash from rs${rsid} : ${hashes.toList().toSet()}")
            return false
        }
    }
}
