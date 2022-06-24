package src

import org.springframework.batch.item.ItemWriter
import org.springframework.data.mongodb.core.MongoTemplate
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

class SSDeprecationWriter implements ItemWriter<SubmittedVariantEntity> {
    MongoTemplate mongoTemplate
    Long accessioningMonotonicInitSs

    SSDeprecationWriter() {

    }

    SSDeprecationWriter(MongoTemplate mongoTemplate, accessioningMonotonicInitSs) {
        this.mongoTemplate = mongoTemplate
        this.accessioningMonotonicInitSs = accessioningMonotonicInitSs
    }

    @Override
    def write(List<? extends SubmittedVariantEntity> svesToDeprecate) {
        def svesToDeprecateInSVE = svesToDeprecate.findAll{sve -> (sve.getAccession() >= accessioningMonotonicInitSs)}
        def svesToDeprecateInDbsnpSVE = svesToDeprecate.findAll{sve -> (sve.getAccession() < accessioningMonotonicInitSs)}
        deprecateVariants(svesToDeprecateInSVE, SubmittedVariantEntity.class)
        deprecateVariants(svesToDeprecateInDbsnpSVE, DbsnpSubmittedVariantEntity.class)
    }

    def deprecateVariants(List<? extends SubmittedVariantEntity> svesToDeprecate, Class sveCollectionToUse) {
        def idsToRemove = svesToDeprecate.collect {sve -> sve.getId()}
        this.mongoTemplate.findAllAndRemove(query(where("_id").in(idsToRemove)), sveCollectionToUse)
        def svoeCollectionToUse = sveCollectionToUse.equals(SubmittedVariantEntity.class) ?
                SubmittedVariantOperationEntity.class : DbsnpSubmittedVariantOperationEntity.class
        writeDeprecationOperation(svesToDeprecate, svoeCollectionToUse)
    }

    def writeDeprecationOperation(List<? extends SubmittedVariantEntity> svesToDeprecate, Class svoeCollectionToUse) {
        def svoesToWrite = svesToDeprecate.collect { sve ->
            SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity()
            svoe.fill(EventType.DEPRECATED, sve.getAccession(), null,
                    "EVA2889: Deprecated since this variant was incorrectly remapped",
                    Collections.singletonList(new SubmittedVariantInactiveEntity(sve)))
            svoe.setId("EVA2889_SS_DEPRECATED_${sve.getId()}")
            return svoe
        }
        this.mongoTemplate.insert(svoesToWrite, svoeCollectionToUse)
    }

}
