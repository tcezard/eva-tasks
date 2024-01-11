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
        def mergedCandidateOperations = []
        is.eachLine{
            String[] sp_line =  it.split("\\t")
            def hash = sp_line.first()
            def rsids = sp_line.tail().collect({it.toLong()}).toArray(Long[]::new)
            mergedCandidateOperations.add(createSingleMergeCandidate(hash, rsids))
            if (mergedCandidateOperations.size() == batchSize){
                totalInsertedCount = insertOperations(mergedCandidateOperations, totalInsertedCount)
                mergedCandidateOperations = [];
            }
        }

        if (mergedCandidateOperations.size() > 0){
            totalInsertedCount = insertOperations(mergedCandidateOperations, totalInsertedCount)
            mergedCandidateOperations = [];
        }
        is.close()
    }

    def insertOperations = {SubmittedVariantOperationEntity[] mergedCandidateOperations, int totalInsertedCount ->
        def insertedCount = devEnv.bulkInsertIgnoreDuplicates(mergedCandidateOperations, SubmittedVariantOperationEntity.class)
        totalInsertedCount = totalInsertedCount + insertedCount
        logger.info("Inserted ${totalInsertedCount} Merge candidate documents")
        return totalInsertedCount
    }

    def createSingleMergeCandidate = { String hash, Long[] rsids ->
        def evaAndDbsnpSveCursors = [sveClass, dbsnpSveClass].collect { collectionClass ->
            prodEnv.mongoTemplate.find(query(where("seq").is(assembly).and("rs").in(rsids)), collectionClass)
        }
        def allSves = evaAndDbsnpSveCursors.each { cursor ->
            cursor.each { sves -> sves}}.flatten()
        def hashSet = allSves.collect { SubmittedVariantEntity sve -> EVAObjectModelUtils.getClusteredVariantHash(sve) }.toSet()
        assertMergeIsValid(allSves, hash)
        // create MERGED_CANDIDATE operation
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

    def assertMergeIsValid = { ArrayList<SubmittedVariantEntity> allSves, String hash ->
        //Check that all the sve provided have the same clustered variant hash which match the one from the input
        def hashSet = allSves.collect { sve -> EVAObjectModelUtils.getClusteredVariantHash(sve) }.toSet()
        assert hashSet.size() == 1
        assert hashSet.head() == hash
    }


}
