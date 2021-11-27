import com.mongodb.MongoCursorNotFoundException
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.0-SNAPSHOT')
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
import org.springframework.data.mongodb.core.query.CriteriaDefinition
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity

import java.time.Instant
import java.time.ZonedDateTime
import java.util.stream.Collectors

import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantEntity

// Verify if any deprecated RS has already been used while carrying out clustering in EVA-2633 and EVA-2634
@Component
@Import(value=[MongoConfiguration.class])
class MainApp implements CommandLineRunner {

    @Autowired
    private MongoTemplate mongoTemplate

    void run(String... args) {

        int numRSCollected = 0
        int batchIndex = 0
        String lastSeenID = null

        while(true) {
            ImmutablePair<List<Document>, String> declusteredRSBatchAndLastSeenID = getNextBatchOfDeprecatedIds(lastSeenID)
            if (declusteredRSBatchAndLastSeenID != null) {
                checkDeprecatedRSUsedByClustering(declusteredRSBatchAndLastSeenID.left)
                numRSCollected += declusteredRSBatchAndLastSeenID.left.size()
                lastSeenID = declusteredRSBatchAndLastSeenID.right
                System.out.println("Checking deprecated RS used by clustering - processed " + numRSCollected + " records...")
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
    void checkDeprecatedRSUsedByClustering(List<Document> declusteredRSToCheck) {
        List<Long> declusteredRSIDs = declusteredRSToCheck*.getLong("accession")
        List<String> assembliesClusteredSoFarIn2021 = Arrays.asList("GCA_000219495.2", "GCA_000372685.2", "GCA_001577835.1", "GCA_000188115.3", "GCA_000146605.4", "GCA_000309985.1", "GCA_001465895.2", "GCA_000181335.4")

        // Assemblies clustered in EVA-2633 and EVA-2634
        CriteriaDefinition assemblyCriteriaInDbsnpCVE = where("asm").in(assembliesClusteredSoFarIn2021)
        Query queryToCheckDeclusteredRSIDsInDbsnpCVECollection =
                query(where("accession").in(declusteredRSIDs)).addCriteria(assemblyCriteriaInDbsnpCVE)
        def declusteredRSInDbsnpCVECollection =  mongoTemplate.find(queryToCheckDeclusteredRSIDsInDbsnpCVECollection, DbsnpClusteredVariantEntity.class)

        // If there were unclustered RS (only possible in EVA collection because all SS in dbSNPSVE collection are already clustered)
        // that were clustered during this exercise or in the past (cannot know for sure due to EVA-2674) identify them
        CriteriaDefinition assemblyCriteriaInSVE = where("seq").in(assembliesClusteredSoFarIn2021)
        List<Long> declusteredRSIDsInDbsnpCVECollection = declusteredRSInDbsnpCVECollection*.getAccession()
        Query queryToCheckDeclusteredRSInEVASVECollection = query(where("rs").in(declusteredRSIDsInDbsnpCVECollection)).addCriteria(assemblyCriteriaInSVE)
        List<SubmittedVariantEntity> ssWithDeclusteredRSIDsInEVASVECollection =
                mongoTemplate.find(queryToCheckDeclusteredRSInEVASVECollection, SubmittedVariantEntity.class)

        ssWithDeclusteredRSIDsInEVASVECollection.forEach{ss -> System.out.println("Error: SS ID " + ss.getAccession() + " has been clustered with declustered RS ID " +  ss.getClusteredVariantAccession() + "!")}
    }
}
