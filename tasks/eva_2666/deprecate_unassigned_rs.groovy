@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.1')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-deprecate', version = '0.6.1')
import com.mongodb.MongoCursorNotFoundException
import org.apache.commons.collections.CollectionUtils
import org.apache.commons.lang3.tuple.ImmutablePair
import org.bson.BsonDocument
import org.bson.BsonInt32
import org.bson.BsonString
import org.bson.BsonValue
import org.bson.Document
import org.bson.conversions.Bson
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.context.annotation.ComponentScan
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import org.springframework.scripting.groovy.GroovyObjectCustomizer
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.deprecate.batch.io.DeprecationWriter

import java.util.stream.Collectors

import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantEntity


@Component
@Import(value=[MongoConfiguration.class])
class DeprecateUnassignedRS implements CommandLineRunner {

    @Autowired
    private MongoTemplate mongoTemplate

    void run(String... args) {

        // RS IDs identified by running check_deprecated_rs_present_in_dbsnpcve.groovy
        List<String> fileLines = new File("eva2666_declustered_rs_ids_to_deprecate.txt").collect {it}
        List<Long> rsIDsToDeprecate = fileLines.collect{Long.parseLong(it)}
        DeprecationWriter deprecationWriter = new DeprecationWriter(mongoTemplate)
        List<Long> rsIDsToDeprecateBatch1 = rsIDsToDeprecate.subList(0, 1000)
        List<Long> rsIDsToDeprecateBatch2 = rsIDsToDeprecate.subList(1000, rsIDsToDeprecate.size())
        deprecationWriter.write(mongoTemplate.find(query(where("accession").in(rsIDsToDeprecateBatch1)), DbsnpClusteredVariantEntity.class))
        deprecationWriter.write(mongoTemplate.find(query(where("accession").in(rsIDsToDeprecateBatch2)), DbsnpClusteredVariantEntity.class))
    }
}