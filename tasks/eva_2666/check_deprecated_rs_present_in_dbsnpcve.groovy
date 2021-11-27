@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.1')
import com.mongodb.MongoCursorNotFoundException
import com.mongodb.client.MongoClient
import org.apache.commons.collections.CollectionUtils
import org.apache.commons.lang3.tuple.ImmutablePair
import org.bson.BsonDocument
import org.bson.BsonInt32
import org.bson.BsonString
import org.bson.BsonValue
import org.bson.Document
import org.bson.conversions.Bson
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity

import java.util.stream.Collectors

import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantEntity


@Import(value=[MongoConfiguration.class])
@Component
class CheckDeprecatedRS implements CommandLineRunner {

    @Autowired
    private MongoTemplate mongoTemplate

    void run(String... args) {

        int numRSCollected = 0
        int batchIndex = 0
        String lastSeenID = null

        while(true) {
            ImmutablePair<List<Document>, String> declusteredRSBatchAndLastSeenID = getNextBatchOfDeprecatedIds(lastSeenID)
            if (declusteredRSBatchAndLastSeenID != null) {
                Collection deprecatedRSAppearingInMainCollection = getDeprecatedRS(mongoTemplate, declusteredRSBatchAndLastSeenID.left)
                deprecatedRSAppearingInMainCollection.forEach{a -> System.out.println("Error: RS ID " + a + " appears in main dbsnpCVE collection but no SS ID is clustered with that RS ID!")}
                numRSCollected += declusteredRSBatchAndLastSeenID.left.size()
                lastSeenID = declusteredRSBatchAndLastSeenID.right
                System.out.println("Checking deprecated RS - processed " + numRSCollected + " records...")
                batchIndex += 1
            } else {
                break
            }
        }
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    ImmutablePair<List<Document>, String> getNextBatchOfDeprecatedIds(String lastSeenID) {
        List<Document> result = new ArrayList<>()

        BsonDocument findQuery = new BsonDocument()
        BsonDocument sortQuery = new BsonDocument("_id", new BsonInt32(1))
        if (Objects.nonNull(lastSeenID)) {
            findQuery = new BsonDocument("_id", new BsonDocument("\$gt", new BsonString(lastSeenID)))
        }
        def declusteredRS = mongoTemplate.getCollection("dbsnpClusteredVariantEntityDeclustered").find(findQuery).sort(sortQuery).limit(1000)
        declusteredRS.iterator().forEachRemaining{a -> result.add(a)}
        if (result.size() > 0) {
            return new ImmutablePair<>(result, result.get(result.size() - 1).getString("_id"))
        }
        return null
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    public Collection getDeprecatedRS(MongoTemplate mongoTemplate, List < Document > declusteredRSToCheck) {
        List<String> delusteredRSHashes = declusteredRSToCheck*.getString("_id")
        List<Long> declusteredRSIDs = declusteredRSToCheck*.getLong("accession")
        Query queryToCheckDeclusteredRSHashesInDbsnpCVECollection = query(where("_id").in(delusteredRSHashes))
        Query queryToCheckDeclusteredRSIDsInDbsnpCVECollection = query(where("accession").in(declusteredRSIDs))
        def declusteredRSInDbsnpCVECollection =  mongoTemplate.find(queryToCheckDeclusteredRSHashesInDbsnpCVECollection, DbsnpClusteredVariantEntity.class)
        declusteredRSInDbsnpCVECollection.addAll(mongoTemplate.find(queryToCheckDeclusteredRSIDsInDbsnpCVECollection, DbsnpClusteredVariantEntity.class))
        List<Long> declusteredRSIDsInDbsnpCVECollection = declusteredRSInDbsnpCVECollection*.getAccession()
        Query queryToCheckDeclusteredRSInDbsnpSVECollection = query(where("rs").in(declusteredRSIDsInDbsnpCVECollection))
        List<Long> declusteredRSIDsInAllSVECollections = mongoTemplate.find(queryToCheckDeclusteredRSInDbsnpSVECollection, DbsnpSubmittedVariantEntity.class)*.getClusteredVariantAccession()
        declusteredRSIDsInAllSVECollections.addAll(mongoTemplate.find(queryToCheckDeclusteredRSInDbsnpSVECollection, SubmittedVariantEntity.class)*.getClusteredVariantAccession())

        // If an RS ID from the declustered collection dbsnpClusteredVariantEntityDeclustered appears in the main
        // dbsnpCVE collection, the only reason must be because there is at least one SS ID clustered with that RS ID
        // Report, if otherwise
        return CollectionUtils.subtract(new HashSet(declusteredRSIDsInDbsnpCVECollection), new HashSet(declusteredRSIDsInAllSVECollections))
    }
}
