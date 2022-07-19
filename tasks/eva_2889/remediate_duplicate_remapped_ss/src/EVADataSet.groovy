package src

import com.mongodb.MongoCursorNotFoundException
import org.apache.commons.lang3.tuple.ImmutablePair
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.data.domain.Sort
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.aggregation.Aggregation
import org.springframework.data.mongodb.core.query.CriteriaDefinition
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

class EVADataSet<T,T2> implements Iterable<T> {
    CriteriaDefinition filterCriteria
    Aggregation aggregation
    MongoTemplate mongoTemplate
    final Class<T> collectionClass
    String pagingAttribute
    final Class<T2> pagingAttributeType
    List<T> currentResult = new ArrayList<>()

    private int numEntriesScanned = 0
    private int batchIndex = 0
    private T2 lastSeenID = null

    private static final Logger logger = LoggerFactory.getLogger(EVADataSet.class)

    EVADataSet() {}

    EVADataSet(CriteriaDefinition filterCriteria, MongoTemplate mongoTemplate, Class<T> collectionClass,
               String pagingAttribute = "_id", Class<T2> pagingAttributeType = String.class) {
        this.filterCriteria = filterCriteria
        this.mongoTemplate = mongoTemplate
        this.collectionClass = collectionClass
        this.pagingAttribute = pagingAttribute
        this.pagingAttributeType = pagingAttributeType
    }

    EVADataSet(Aggregation aggregation, MongoTemplate mongoTemplate, Class<T> collectionClass,
               String pagingAttribute = "_id") {
        this.aggregation = aggregation
        this.mongoTemplate = mongoTemplate
        this.collectionClass = collectionClass
        this.pagingAttribute = pagingAttribute
    }

    @Override
    public Iterator<List<T>> iterator() {
        return new Iterator<List<T>>() {
            @Override
            boolean hasNext() {
                getNextBatch()
                return !(Objects.isNull(this.currentResult))
            }

            @Override
            List<T> next() {
                return this.currentResult
            }
        }
    }

    public void reset() {
        this.lastSeenID = null
        this.batchIndex = 0
        this.numEntriesScanned = 0
    }

    public void getNextBatch() {
        ImmutablePair<List<T>, T2> resultPair = Objects.isNull(this.aggregation)? getNextBatchWithFind(): null
        if (resultPair != null) {
            this.numEntriesScanned += resultPair.left.size()
            this.lastSeenID = resultPair.right
            this.batchIndex += 1
            logger.info("${this.numEntriesScanned} entries of ${collectionClass.getSimpleName()} collection scanned so far...")
            this.currentResult = resultPair.left
            return
        }
        this.currentResult = null
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    private ImmutablePair<List<T>, T2> getNextBatchWithFind() {
        Query queryToGetNextBatch = new Query()
        if (Objects.nonNull(this.filterCriteria)) {
            queryToGetNextBatch.addCriteria(this.filterCriteria)
        }
        if (Objects.nonNull(this.lastSeenID)) {
            queryToGetNextBatch.addCriteria(where(pagingAttribute).gt(this.lastSeenID))
        }
        queryToGetNextBatch = queryToGetNextBatch.with(Sort.by(Sort.Direction.ASC, pagingAttribute)).limit(1000)
        def result = this.mongoTemplate.find(queryToGetNextBatch, this.collectionClass)
        if (result.size() > 0) {
            return new ImmutablePair<>(result, result.get(result.size() - 1).getProperties()
                    .get(this.pagingAttribute.equals("_id")? "id":this.pagingAttribute))
        }
        return null
    }

    public void copyTo(EVADatabaseEnvironment targetDBEnv) {
        this.each {List<T> dbEntities ->
            def entitiesToInsert = dbEntities.groupBy {dbEntity -> dbEntity.getId()};
            def alreadyExistingInTarget = targetDBEnv.mongoTemplate.find(query(where("_id").in(entitiesToInsert.keySet())),
                    this.collectionClass).collect{dbEntity -> dbEntity.getId()}.toSet();
            entitiesToInsert = entitiesToInsert.findAll {k, v -> !alreadyExistingInTarget.contains(k)};
            targetDBEnv.mongoTemplate.insert(entitiesToInsert.values().flatten(), this.collectionClass)
        }
    }
}
