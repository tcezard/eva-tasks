package src

@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-clustering', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-deprecate', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'variation-commons-batch', version = '0.8.1')

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

    void run(String... args) {

        def prodRSFromSVOE =  new EVADataSet(where("_id").regex("SS_DEPRECATED.*EVA2889.*"), this.mongoTemplate,
                SubmittedVariantOperationEntity.class);
        def prodRSDeprecatedCvoe = new EVADataSet(where("_id").regex("RS_DEPRECATED.*EVA2889.*"), this.mongoTemplate,
                ClusteredVariantOperationEntity.class, "accession", Long.class);
        def prodRSDeprecatedDbsnpCvoe = new EVADataSet(where("_id").regex("RS_DEPRECATED.*EVA2889.*"), this.mongoTemplate,
                DbsnpClusteredVariantOperationEntity.class, "accession", Long.class);
        def checkPresentInSubmittedOperations = { cvoes ->
            def deprecatedRsIDs = cvoes.collect{it.getInactiveObjects().get(0).getAccession()};
            def rsIDsFromSvoeAndDbsnpSvoe = (this.mongoTemplate.find(query(where("inactiveObjects.rs").in(deprecatedRsIDs)),
                    SubmittedVariantOperationEntity.class).collect{it.getInactiveObjects().get(0).getClusteredVariantAccession()} +
                    this.mongoTemplate.find(query(where("inactiveObjects.rs").in(deprecatedRsIDs)),
                            DbsnpSubmittedVariantOperationEntity.class).collect{it.getInactiveObjects().get(0).getClusteredVariantAccession()})
                    .unique();
            println("Checking ${deprecatedRsIDs.size()} deprecated RS IDs in SVOE and dbsnpSVOE...");
            println("Found ${rsIDsFromSvoeAndDbsnpSvoe.size()} deprecated RS IDs in SVOE and dbsnpSVOE...");
            (deprecatedRsIDs - rsIDsFromSvoeAndDbsnpSvoe).each {println("ERROR: Found deprecated RS ID ${it} not associated with deprecated SS...")};
        }
        def checkPresentInSubmittedEntities = { cvoes ->
            cvoes.groupBy{it.getInactiveObjects().get(0).getAssemblyAccession()}.each {String assembly, cvoesGrouped ->
                List<Long> deprecatedRsIDs = cvoesGrouped.collect { it.getInactiveObjects().get(0).getAccession() };
                def rsIDsFromSveAndDbsnpSve =
                        this.submittedVariantAccessioningService.getByClusteredVariantAccessionIn(deprecatedRsIDs)
                                .findAll{it.getData().getReferenceSequenceAccession().equals(assembly)}
                                .collect{it.getData().getClusteredVariantAccession()};
                println("Checking ${deprecatedRsIDs.size()} deprecated RS IDs in SVE and dbsnpSVE...");
                println("Found ${rsIDsFromSveAndDbsnpSve.size()} deprecated RS IDs in SVE and dbsnpSVE...");
                rsIDsFromSveAndDbsnpSve.each { println("ERROR: Found deprecated RS ID ${it} associated with current SS...") };
            }
        }
        // Ensure that all RS recorded as deprecated is present in SVOE or dbsnpSVOE - i.e., they were due to deprecated SS
        prodRSDeprecatedCvoe.each {checkPresentInSubmittedOperations(it)};
        prodRSDeprecatedDbsnpCvoe.each {checkPresentInSubmittedOperations(it)};

        // Ensure that none of the SS have any of the deprecated RS assigned to them
        prodRSDeprecatedCvoe.reset(); prodRSDeprecatedDbsnpCvoe.reset();
        prodRSDeprecatedCvoe.each {checkPresentInSubmittedEntities(it)};
        prodRSDeprecatedDbsnpCvoe.each {checkPresentInSubmittedEntities(it)};
    }
}
