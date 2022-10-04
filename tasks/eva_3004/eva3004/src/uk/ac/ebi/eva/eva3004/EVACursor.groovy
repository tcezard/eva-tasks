package uk.ac.ebi.eva.eva3004

import com.mongodb.MongoCursorNotFoundException
import com.mongodb.client.MongoCursor
import org.bson.Document
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.CriteriaDefinition
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable

// Implement a batching retryable cursor on top of MongoCursor
class EVACursor<T> implements Iterable<T> {
    CriteriaDefinition filterCriteria
    MongoTemplate mongoTemplate
    final Class<T> collectionClass
    final int pageSize
    MongoCursor<Document> resultIterator

    EVACursor() {}

    EVACursor(CriteriaDefinition filterCriteria, MongoTemplate mongoTemplate, Class<T> collectionClass, int pageSize = 1000) {
        this.filterCriteria = filterCriteria
        this.mongoTemplate = mongoTemplate
        this.collectionClass = collectionClass
        this.pageSize = pageSize

        this.resultIterator = this.mongoTemplate.getCollection(this.mongoTemplate
                .getCollectionName(this.collectionClass))
                .find(this.filterCriteria.criteriaObject).noCursorTimeout(true).batchSize(pageSize).iterator()
    }

    @Override
    public Iterator<List<T>> iterator() {
        return new Iterator<List<T>>() {
            @Override
            @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
            boolean hasNext() {
                this.resultIterator.hasNext()
            }

            @Override
            List<T> next() {
                List<T> result = new ArrayList<>()
                result.add(this.mongoTemplate.converter.read(this.collectionClass, this.resultIterator.next()))
                for (i in 0..this.pageSize-2) {
                    if (this.resultIterator.hasNext()) {
                        result.add(this.mongoTemplate.converter.read(this.collectionClass, this.resultIterator.next()))
                    } else {
                        break
                    }
                }
                return result
            }
        }
    }
}
