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

class CreateMergeCandidates {
    static def logger = LoggerFactory.getLogger(CreateMergeCandidates.class)

    String assembly
    String mergeCandidatesFile
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv
    int batchSize

    CreateMergeCandidates() {}

    CreateMergeCandidates(String assembly, String mergeCandidatesFile, EVADatabaseEnvironment prodEnv,
                          EVADatabaseEnvironment devEnv, int batchSize) {
        this.assembly = assembly
        this.mergeCandidatesFile = mergeCandidatesFile
        this.prodEnv = prodEnv
        this.devEnv = devEnv
        this.batchSize = batchSize
    }


    def process = {
        def is = new File(mergeCandidatesFile).newInputStream()
        logger.info("Create Merge candidate documents in batches of ${batchSize}")
        int totalInsertedCount = 0
        int lineRead = 0
        def mergedCandidateOperations = []
        ArrayList<String> batchedLines = new ArrayList<String>()
        is.eachLine{
            lineRead = lineRead + 1
            batchedLines.add(it)
            if (batchedLines.size() == batchSize){
                def mergeCandidateOperations = createBatchMergeCandidates(batchedLines)
                totalInsertedCount = insertOperations(mergeCandidateOperations, totalInsertedCount, lineRead)
                batchedLines = new ArrayList()
            }
        }
        if (batchedLines.size() > 0){
            def mergeCandidateOperations = createBatchMergeCandidates(batchedLines)
            totalInsertedCount = insertOperations(mergeCandidateOperations, totalInsertedCount, lineRead)
        }
        is.close()

    }

    def createBatchMergeCandidates = { List<String> hashAndRsids ->
        Map<String, Long[]> hashToRsids = hashAndRsids.collectEntries {
            String[] sp_line =  it.split("\\t")
            def hash = sp_line.first()
            def rsids = sp_line.tail().collect({it.toLong()})
            return [ (hash): rsids]
        }
        def rsids = hashToRsids.values().toList().flatten()
        def allSves = [sveClass, dbsnpSveClass]
                .collectMany( collectionClass ->
            prodEnv.mongoTemplate.find(query(where("seq").is(assembly).and("rs").in(rsids)), collectionClass)
                ).flatten()
        Map<String, SubmittedVariantEntity[]> hashToSves = new HashMap<>()
        allSves.each {SubmittedVariantEntity sve ->
            def hash = EVAObjectModelUtils.getClusteredVariantHash(sve)
            if (! hashToSves.containsKey(hash) ){
                hashToSves.put(hash, new ArrayList<SubmittedVariantEntity>())
            }
            hashToSves[hash].add(sve)
        }
        def mergeCandidateOperations= hashToSves.keySet().collect {String hash ->
            allSves = hashToSves.get(hash)
            rsids = allSves.collect(sve -> sve.clusteredVariantAccession).toSet().toList()
            if (! isMergeIsValid(allSves, hash, rsids)) {
                return null
            }
            // create MERGE_CANDIDATE operation
            def submittedVariantInactiveEntity = allSves.collect {new SubmittedVariantInactiveEntity(it)}
            SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity()
            Long rsAccession = ((SubmittedVariantEntity)allSves[0]).clusteredVariantAccession
            svoe.fill(
                    EventType.RS_MERGE_CANDIDATES, rsAccession,
                    "RS mismatch with " + rsAccession, submittedVariantInactiveEntity
            )
            svoe.setId(String.format("%s_%s_%s", RSMergeAndSplitCandidatesReaderConfiguration.MERGE_CANDIDATE_ID_PREFIX,
                                     assembly, hash))
            return svoe
        }
        // Remove null for invalid MERGE candidates
        mergeCandidateOperations.removeAll([null])
        return mergeCandidateOperations
    }

    def insertOperations = {ArrayList mergedCandidateOperations, int totalInsertedCount, int lineRead ->
        def insertedCount = devEnv.bulkInsertIgnoreDuplicates(mergedCandidateOperations, SubmittedVariantOperationEntity.class)
        totalInsertedCount = totalInsertedCount + insertedCount
        logger.info("Read ${lineRead}, Inserted ${totalInsertedCount} Merge candidate documents")
        return totalInsertedCount
    }

    def isMergeIsValid = { ArrayList<SubmittedVariantEntity> allSves, String hash, ArrayList<Long> rsids->
        //Check that all the sve provided have the same clustered variant hash which match the one from the input
        def hashSet = allSves.collect { sve -> EVAObjectModelUtils.getClusteredVariantHash(sve) }.toSet()
        if (hashSet.size() == 1 && hashSet.head() == hash){
            return true
        }else if (hashSet.size() > 1) {
            logger.warn("More than one Hash in associated with the group of rsids ${rsids}")
            return false
        }else{
            logger.warn("Hash from file ${hash} is different from Hash from the one in the one in database ${hashSet.head()}")
            return false
        }
    }


}
