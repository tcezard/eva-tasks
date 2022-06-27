@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')

import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.boot.CommandLineRunner
import org.springframework.context.annotation.Import
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

@Component
@Import(value=[MongoConfiguration.class])
class QCRemediation implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(QCRemediation.class)

    @Autowired
    private MongoTemplate mongoTemplate

    @Value("${assembly}")
    private String assembly

    @Value("${unremappedVcf}")
    private String unremappedVcf
    

    void run(String... args) {
        Query getOperationsQuery = query(where("_id").regex("^SS_SPLIT_.+").and("inactiveObjects.seq").is(assembly))
        List<? extends SubmittedVariantOperationEntity> operations = this.mongoTemplate.find(getOperationsQuery, DbsnpSubmittedVariantOperationEntity.class)
        List<Long> accessions = operations.collect{ it.getSplitInto() }

        File unmappedFile = new File(unremappedVcf)

        // These accessions should be in either unremapped file, or EVA collection, but not in dbsnp
        for (Long accession : accessions) {
            if (!unmappedFile.text.contains('ss' + accession)) {
                Query findRemappedByAccession = query(where("accession").is(accession).and("remappedFrom").exists(true))
                List<DbsnpSubmittedVariantEntity> dbsnpEntities = this.mongoTemplate.find(findRemappedByAccession, DbsnpSubmittedVariantEntity.class)
                List<SubmittedVariantEntity> evaEntities = this.mongoTemplate.find(findRemappedByAccession, SubmittedVariantEntity.class)
                if (dbsnpEntities.size() > 0) {
                    logger.error("Found " + dbsnpEntities.size() + " dbSNP entries with accession " + accession)
                }
                if (evaEntities.size() == 0) {
                    logger.error("Found no EVA entries with accession " + accession)
                }
            }
        }
    }

}
