package uk.ac.ebi.eva.eva3004

import com.mongodb.MongoClient
import com.mongodb.MongoClientOptions
import com.mongodb.ReadConcern
import com.mongodb.ReadPreference
import com.mongodb.WriteConcern
import org.bson.Document
import org.springframework.boot.autoconfigure.mongo.MongoClientFactory
import org.springframework.boot.autoconfigure.mongo.MongoProperties
import org.springframework.context.annotation.AnnotationConfigApplicationContext
import org.springframework.core.convert.converter.Converter
import org.springframework.core.env.PropertiesPropertySource
import org.springframework.core.env.StandardEnvironment
import org.springframework.dao.DuplicateKeyException
import org.springframework.data.mapping.model.SimpleTypeHolder
import org.springframework.data.mongodb.core.BulkOperations
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.SimpleMongoDbFactory
import org.springframework.data.mongodb.core.convert.DefaultDbRefResolver
import org.springframework.data.mongodb.core.convert.DefaultMongoTypeMapper
import org.springframework.data.mongodb.core.convert.MappingMongoConverter
import org.springframework.data.mongodb.core.mapping.MongoMappingContext
import org.springframework.data.mongodb.core.mapping.MongoSimpleTypes
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.model.ClusteredVariant
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.service.nonhuman.ClusteredVariantAccessioningService
import uk.ac.ebi.eva.accession.core.service.nonhuman.SubmittedVariantAccessioningService
import uk.ac.ebi.eva.accession.core.summary.ClusteredVariantSummaryFunction
import uk.ac.ebi.eva.commons.core.models.VariantClassifier

import java.time.LocalDateTime
import java.time.ZoneId

import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where

// This class is to create an environment for EVA databases from a given properties file
// We do this instead of Spring Boot autowiring because Spring boot
// doesn't support multiple application contexts to hold properties
// for multiple environments (ex: prod and dev) at the same time
// Most of the code below is borrowed from MongoConfiguration in the eva-accession repository
public class EVADatabaseEnvironment {
    MongoClient mongoClient
    MongoTemplate mongoTemplate
    SubmittedVariantAccessioningService submittedVariantAccessioningService
    ClusteredVariantAccessioningService clusteredVariantAccessioningService
    AnnotationConfigApplicationContext springApplicationContext

    static def sveClass = SubmittedVariantEntity.class
    static def dbsnpSveClass = DbsnpSubmittedVariantEntity.class
    static def svoeClass = SubmittedVariantOperationEntity.class
    static def dbsnpSvoeClass = DbsnpSubmittedVariantOperationEntity.class
    static def cveClass = ClusteredVariantEntity.class
    static def dbsnpCveClass = DbsnpClusteredVariantEntity.class
    static def cvoeClass = ClusteredVariantOperationEntity.class
    static def dbsnpCvoeClass = DbsnpClusteredVariantOperationEntity.class

    EVADatabaseEnvironment(MongoClient mongoClient, MongoTemplate mongoTemplate,
                           SubmittedVariantAccessioningService submittedVariantAccessioningService,
                           ClusteredVariantAccessioningService clusteredVariantAccessioningService,
                           AnnotationConfigApplicationContext springApplicationContext) {
        this.mongoClient = mongoClient
        this.mongoTemplate = mongoTemplate
        this.submittedVariantAccessioningService = submittedVariantAccessioningService
        this.clusteredVariantAccessioningService = clusteredVariantAccessioningService
        this.springApplicationContext = springApplicationContext
    }

    EVADatabaseEnvironment(MongoClient mongoClient, MongoTemplate mongoTemplate) {
        this.mongoClient = mongoClient
        this.mongoTemplate = mongoTemplate
    }

    static EVADatabaseEnvironment createFromSpringContext(String propertiesFile, Class springApplicationClass) {
        AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext()
        def appProps = new Properties()
        appProps.load(new FileInputStream(new File(propertiesFile)))
        context.getEnvironment().getPropertySources().addLast(new PropertiesPropertySource("main", appProps))
        context.register(springApplicationClass)
        context.refresh()

        def mc = context.getBean(MongoClient.class)
        def mt = context.getBean(MongoTemplate.class)
        def sva = context.getBean(SubmittedVariantAccessioningService.class)
        def cva = context.getBean(ClusteredVariantAccessioningService.class)
        return new EVADatabaseEnvironment(mc, mt, sva, cva, context)
    }

    static EVADatabaseEnvironment parseFrom(String propertiesFile) {
        def properties = new Properties()
        properties.load(new FileInputStream(new File(propertiesFile)))
        MongoProperties mongoProperties = getMongoPropertiesForEnv(properties)
        MongoClient mongoClient = getMongoClientForEnv(properties, mongoProperties)
        MongoTemplate mongoTemplate = getMongoTemplateForEnv(mongoClient, mongoProperties)

        return new EVADatabaseEnvironment(mongoClient, mongoTemplate)
    }

    private static MongoClient getMongoClientForEnv(Properties properties, MongoProperties mongoProperties) {
        def readPreference = ReadPreference.valueOf(properties.getProperty("mongodb.read-preference"))
        def mongoClientOptions = new MongoClientOptions.Builder().readPreference(readPreference).writeConcern(WriteConcern.MAJORITY).readConcern(ReadConcern.MAJORITY)

        def environment = new StandardEnvironment()
        def mongoClient = new MongoClientFactory(mongoProperties, environment).createMongoClient(mongoClientOptions.build())
        mongoClient
    }

    private static MongoProperties getMongoPropertiesForEnv(Properties properties) {
        def mongoProperties = new MongoProperties()
        mongoProperties.setDatabase(properties.getProperty("spring.data.mongodb.database"))
        mongoProperties.setHost(properties.getProperty("spring.data.mongodb.host"))
        mongoProperties.setPort(Integer.parseInt(properties.getProperty("spring.data.mongodb.port")))
        mongoProperties.setUsername(properties.getProperty("spring.data.mongodb.username"))
        mongoProperties.setPassword(properties.getProperty("spring.data.mongodb.password").toCharArray())
        mongoProperties.setAuthenticationDatabase(properties.getProperty("spring.data.mongodb.authentication-database"))
        return mongoProperties
    }

    private static MongoTemplate getMongoTemplateForEnv(MongoClient mongoClient, MongoProperties mongoProperties) {
        def mongoDbFactory = new SimpleMongoDbFactory(mongoClient, mongoProperties.getDatabase())
        def mappingContext = new MongoMappingContext()
        SimpleTypeHolder simpleTypeHolder = new SimpleTypeHolder(new HashSet<>(Arrays.asList(
                Date.class,
                LocalDateTime.class
        )), MongoSimpleTypes.HOLDER)
        mappingContext.setSimpleTypeHolder(simpleTypeHolder)
        def mappingMongoConverter = new MappingMongoConverter(new DefaultDbRefResolver(mongoDbFactory), mappingContext)
        mappingMongoConverter.setTypeMapper(new DefaultMongoTypeMapper(null))
        List<Converter<?, ?>> converterList = new ArrayList<Converter<?, ?>>()
        converterList.add(new MongoLocalDateTimeFromStringConverter())
        converterList.add(new MongoDateTimeFromStringConverter())
        mappingMongoConverter.setMapKeyDotReplacement("#")
        mappingMongoConverter.afterPropertiesSet()
        return new MongoTemplate(mongoDbFactory, mappingMongoConverter)
    }

    private static final class MongoLocalDateTimeFromStringConverter implements Converter<String, LocalDateTime> {
        @Override
        LocalDateTime convert(String source) {
            return source == null ? null : LocalDateTime.parse(source)
        }
    }

    private static final class MongoDateTimeFromStringConverter implements Converter<String, Date> {
        @Override
        Date convert(String source) {
            return source == null ? null :
                    Date.from(LocalDateTime.parse(source).toLocalDate()
                            .atStartOfDay(ZoneId.systemDefault()).toInstant())
        }
    }

    static def shardCollection = {collectionName, dbEnvironment ->
        def fqdbName = dbEnvironment.mongoTemplate.getDb().name + "." + collectionName
        def shardCommandDocument = new Document("shardCollection", fqdbName)
        shardCommandDocument.put("key", new Document("_id", 1))
        dbEnvironment.mongoClient.getDatabase("admin").runCommand(shardCommandDocument)
    }

    void restoreCollectionsToOriginalState() {
        // Restore collections to their original state
        this.mongoTemplate.getCollectionNames().each { collectionName ->
            def before_remediation_suffix = "_before_remediation"
            def targetCollectionName = collectionName.replace(before_remediation_suffix, "")
            if (collectionName.endsWith(before_remediation_suffix)) {
                this.mongoTemplate.getCollection(targetCollectionName).drop()
                this.mongoTemplate.getCollection(collectionName).aggregate(
                        Collections.singletonList(new Document("\$out",
                                targetCollectionName))).allowDiskUse(true).size()
                shardCollection(targetCollectionName, this)
            }
        }
    }

    private void backup(String backupSuffix) {
        this.mongoTemplate.getCollectionNames().each {collectionName ->
            def before_remediation_suffix = "_before_remediation"
            def after_remediation_suffix = "_after_remediation"
            def targetCollectionName = collectionName + backupSuffix
            if (!collectionName.endsWith(after_remediation_suffix) && !collectionName.endsWith(before_remediation_suffix)) {
                this.mongoTemplate.getCollection(collectionName).aggregate(
                        Collections.singletonList(new Document("\$out",
                                targetCollectionName))).allowDiskUse(true).size()
            }
        }
    }

    public void backupCollectionsAfterRemediation() {
        backup("_after_remediation")
    }

    public void backupCollectionsBeforeRemediation() {
        backup("_before_remediation")
    }

    public <T> void bulkInsertIgnoreDuplicates(List<T> recordsToInsert, Class<T> collectionClass, String collectionName = null) {
        def recordsInsertedResult = null
        if (recordsToInsert.size() > 0) {
            BulkOperations ops
            if (Objects.isNull(collectionName)) {
                ops = this.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, collectionClass)
            } else {
                ops = this.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, collectionClass, collectionName)
            }
            ops.insert(recordsToInsert)
            try {
                ops.execute()
            }
            catch(DuplicateKeyException ignored) {

            }
        }
    }

    static def getClusteredVariantHash = {submittedVariant ->
        ClusteredVariant clusteredVariant = toClusteredVariant(submittedVariant)
        def clusteredHashingFunction = (new ClusteredVariantSummaryFunction()).andThen(new SHA1HashingFunction())
        return clusteredHashingFunction.apply(clusteredVariant)
    }

    static def toClusteredVariant = {submittedVariant ->
        return new ClusteredVariant(submittedVariant.getReferenceSequenceAccession(), submittedVariant.getTaxonomyAccession(), submittedVariant.getContig(), submittedVariant.getStart(), VariantClassifier.getVariantClassification(submittedVariant.getReferenceAllele(), submittedVariant.getAlternateAllele()), submittedVariant.isValidated(), submittedVariant.getCreatedDate())
    }

    static def toClusteredVariantEntity = {submittedVariantEntity ->
        return new ClusteredVariantEntity(submittedVariantEntity.getClusteredVariantAccession(),
                getClusteredVariantHash(submittedVariantEntity), toClusteredVariant(submittedVariantEntity))
    }
}
