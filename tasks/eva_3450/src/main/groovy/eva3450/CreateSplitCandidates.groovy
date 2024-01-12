package eva3450

import com.google.common.collect.Iterables
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

        is.eachLine{
            lineRead = lineRead + 1
            String[] sp_line = it.split("\\t")
            def rsid = sp_line.first().toLong()
            def hashes = sp_line.tail()
            splitCandidateOperations.add(createSingleSplitCandidate(rsid, hashes))
            if (splitCandidateOperations.size() == batchSize){
                totalInsertedCount = insertOperations(splitCandidateOperations, totalInsertedCount, lineRead)
                splitCandidateOperations = [];
            }
        }
        if (splitCandidateOperations.size() > 0){
            totalInsertedCount = insertOperations(splitCandidateOperations, totalInsertedCount, lineRead)
            splitCandidateOperations = [];
        }
        is.close()

    }
    def insertOperations = {ArrayList splitCandidateOperations, int totalInsertedCount, int lineRead ->
        def insertedCount = devEnv.bulkInsertIgnoreDuplicates(splitCandidateOperations, SubmittedVariantOperationEntity.class)
        totalInsertedCount = totalInsertedCount + insertedCount
        logger.info("Read ${lineRead}, Inserted ${totalInsertedCount} Split candidate documents")
        return totalInsertedCount
    }

    def createSingleSplitCandidate = { Long rsid, String[] hashes ->
        def evaAndDbsnpSveCursors = [sveClass, dbsnpSveClass].collect { collectionClass ->
            prodEnv.mongoTemplate.find(query(where("seq").is(assembly).and("rs").is(rsid)), collectionClass)
        }
        def allSves = evaAndDbsnpSveCursors.each { cursor ->
            cursor.each { sves -> sves}}.flatten()
        assertSplitIsValid(allSves, hashes)
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

    def assertSplitIsValid = { ArrayList<SubmittedVariantEntity> allSves, String[] hashes ->
        //Check that all the sve provided have the same clustered variant hash which match the one from the input
        def hashSet = allSves.collect { sve -> EVAObjectModelUtils.getClusteredVariantHash(sve) }.toSet()
        assert hashSet == hashes.toList().toSet()
    }

}
