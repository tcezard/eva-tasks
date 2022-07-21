package src

@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-clustering', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-deprecate', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'variation-commons-batch', version = '0.8.1')

import com.mongodb.ReadPreference
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.boot.CommandLineRunner
import org.springframework.context.annotation.Import
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.accession.clustering.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.io.SSSplitWriterConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.listeners.ListenersConfiguration
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.ClusteredVariantAccessioningConfiguration
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.SubmittedVariantAccessioningConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.*;
import uk.ac.ebi.eva.accession.core.model.eva.*
import uk.ac.ebi.eva.accession.core.service.nonhuman.ClusteredVariantAccessioningService
import uk.ac.ebi.eva.accession.core.service.nonhuman.SubmittedVariantAccessioningService
import uk.ac.ebi.eva.accession.deprecate.configuration.batch.io.SubmittedVariantDeprecationWriterConfiguration
import uk.ac.ebi.eva.metrics.configuration.MetricConfiguration;

import static org.springframework.data.mongodb.core.query.Query.query;
import static org.springframework.data.mongodb.core.query.Criteria.where;

@Component
@Import(value=[SSSplitWriterConfiguration.class,
        InputParametersConfiguration.class, MetricConfiguration.class, ListenersConfiguration.class,
        SubmittedVariantAccessioningConfiguration.class, ClusteredVariantAccessioningConfiguration.class,
        SubmittedVariantDeprecationWriterConfiguration.class])
class QCRemediatedSS implements CommandLineRunner {

    @Autowired
    private MongoTemplate mongoTemplate

    @Autowired
    private SubmittedVariantAccessioningService submittedVariantAccessioningService

    @Autowired
    private ClusteredVariantAccessioningService clusteredVariantAccessioningService

    void run(String... args) {
        this.mongoTemplate.setReadPreference(ReadPreference.secondaryPreferred())
        def ssDeprecatedOpIDPattern = "SS_DEPRECATED.*EVA2889.*"
        def svoeDeprecationOps =  new EVADataSet(where("_id").regex(ssDeprecatedOpIDPattern), this.mongoTemplate,
                SubmittedVariantOperationEntity.class, "accession", Long.class);
        def dbsnpSvoeDeprecationOps =  new EVADataSet(where("_id").regex(ssDeprecatedOpIDPattern), this.mongoTemplate,
                DbsnpSubmittedVariantOperationEntity.class, "accession", Long.class);


        def rsDeprecatedOpIDPattern = "RS_DEPRECATED.*EVA2889.*"
        def cvoeDeprecationOps = new EVADataSet(where("_id").regex(rsDeprecatedOpIDPattern), this.mongoTemplate,
                ClusteredVariantOperationEntity.class, "accession", Long.class);
        def dbsnpCvoeDeprecationOps = new EVADataSet(where("_id").regex(rsDeprecatedOpIDPattern), this.mongoTemplate,
                DbsnpClusteredVariantOperationEntity.class, "accession", Long.class);

        def getSVEsWithRSIDsInAssembly = {assembly, rsIDs ->
            [SubmittedVariantEntity.class, DbsnpSubmittedVariantEntity.class].collect {collection ->
                this.mongoTemplate.find(query(where("seq").is(assembly).and("rs").in(rsIDs)),collection)
            }.flatten()
        }
        def checkPresentInSubmittedOperations = { cvoes ->
            def deprecatedRsIDs = cvoes.collect{it.getInactiveObjects().get(0).getAccession()}.unique();
            def rsIDsFromSvoeAndDbsnpSvoe =
                    [SubmittedVariantOperationEntity.class, DbsnpSubmittedVariantOperationEntity.class].collect{collection ->
                        this.mongoTemplate.find(query(where("_id").regex(ssDeprecatedOpIDPattern)
                                .and("inactiveObjects.rs").in(deprecatedRsIDs)), collection)
                                .collect{it.getInactiveObjects().get(0).getClusteredVariantAccession()}
                    }.flatten().unique();
            println("Checking ${deprecatedRsIDs.size()} deprecated RS IDs in SVOE and dbsnpSVOE...");
            println("Found ${rsIDsFromSvoeAndDbsnpSvoe.size()} deprecated RS IDs in SVOE and dbsnpSVOE...");
            (deprecatedRsIDs - rsIDsFromSvoeAndDbsnpSvoe).each {println("ERROR: Found deprecated RS ID ${it} not associated with deprecated SS!!")};
        }
        def checkRSNotAssigned = { cvoes ->
            cvoes.groupBy{it.getInactiveObjects().get(0).getAssemblyAccession()}.each {String assembly, cvoesGrouped ->
                List<Long> deprecatedRsIDs = cvoesGrouped.collect { it.getInactiveObjects().get(0).getAccession() }.unique();
                println("Checking ${deprecatedRsIDs.size()} deprecated RS IDs in SVE and dbsnpSVE...");
                def rsIDsFromSveAndDbsnpSve = getSVEsWithRSIDsInAssembly(assembly, deprecatedRsIDs)
                        .collect{it.getData().getClusteredVariantAccession()}.unique();
                println("Found ${rsIDsFromSveAndDbsnpSve.size()} deprecated RS IDs in SVE and dbsnpSVE...");
                rsIDsFromSveAndDbsnpSve.each { println("ERROR: Found deprecated RS ID ${it} associated with current SS!!") };
            }
        }
        def checkSSNotPresentInSVE = {svoes ->
            svoes.groupBy{it.getInactiveObjects().get(0).getReferenceSequenceAccession()}.each {assembly, svoesGrouped ->
                def deprecatedSVEHashes = svoesGrouped.collect {svoe ->
                    svoe.getInactiveObjects().get(0).getHashedMessage()};
                println("Checking ${deprecatedSVEHashes.size()} deprecated SS in SVE and dbsnpSVE...")
                [SubmittedVariantEntity.class, DbsnpSubmittedVariantEntity.class].each {collection ->
                    this.mongoTemplate.find(query(where("_id").in(deprecatedSVEHashes)), collection).each{
                        println("ERROR: Found deprecated SS with hash ${it.getHashedMessage()} in ${collection.getSimpleName()}!!")
                    }
                }
            }
        }
        def checkNonDeprecatedRSNotOrphaned = {svoes ->
            svoes.groupBy{it.getInactiveObjects().get(0).getReferenceSequenceAccession()}.each {assembly, svoesGrouped ->
                def associatedRSIDs = svoesGrouped.collect {svoe ->
                    svoe.getInactiveObjects().get(0).getClusteredVariantAccession()}.unique();
                def deprecatedRSIDs =
                        [ClusteredVariantOperationEntity.class, DbsnpClusteredVariantOperationEntity.class].collect {collection ->
                            return this.mongoTemplate.find(query(where("_id").regex(rsDeprecatedOpIDPattern)
                                    .and("accession").in(associatedRSIDs)), collection).collect{it.getAccession()}
                        }.flatten().unique();
                def nonDeprecatedRSIDs = associatedRSIDs - deprecatedRSIDs;
                def nonDeprecatedRSIDsInSVE = getSVEsWithRSIDsInAssembly(assembly, nonDeprecatedRSIDs)
                        .groupBy {it.getClusteredVariantAccession()}
                (nonDeprecatedRSIDs - nonDeprecatedRSIDsInSVE.keySet()).each {
                    println("ERROR: RS ID ${it} should have been deprecated!!")
                }
            }
        }

        [cvoeDeprecationOps, dbsnpCvoeDeprecationOps].each {ops -> ops.each{
            // Ensure that all RS recorded as deprecated is present in SVOE or dbsnpSVOE - i.e., they were due to deprecated SS
            checkPresentInSubmittedOperations(it);
            // Ensure that none of the SS have any of the deprecated RS assigned to them
            checkRSNotAssigned(it);
        }};

        // Ensure no deprecated SS recorded in dbsnpSVOE and SVOE appear in the Submitted collections
        [svoeDeprecationOps, dbsnpSvoeDeprecationOps].each {ops -> ops.each {
            checkSSNotPresentInSVE(it);
            checkNonDeprecatedRSNotOrphaned(it);
        }};
    }
}
