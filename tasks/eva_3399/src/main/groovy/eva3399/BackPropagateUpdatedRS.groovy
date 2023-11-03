package eva3399

import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.eva.accession.clustering.batch.io.BackPropagatedRSWriter
import uk.ac.ebi.eva.accession.clustering.configuration.BeanNames
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

class BackPropagateUpdatedRS {
    String remappedAssembly
    String sourceAssembly
    EVADatabaseEnvironment backPropEnv

    BackPropagateUpdatedRS(String remappedAssembly, String sourceAssembly, EVADatabaseEnvironment backPropEnv) {
        this.remappedAssembly = remappedAssembly
        this.sourceAssembly = sourceAssembly
        this.backPropEnv = backPropEnv
    }

    void backPropagate() {
        def backPropagatedRSWriter =
                this.backPropEnv.springApplicationContext.getBean(BeanNames.BACK_PROPAGATED_RS_WRITER,
                        BackPropagatedRSWriter.class)
        def opsOfImpactedSS = [svoeClass, dbsnpSvoeClass].collect{
            [new RetryableBatchingCursor<>(
                where("_id").regex("^EVA3399_UPD_LOCUS_${this.remappedAssembly}_.*").and(
                        "inactiveObjects.remappedFrom").is(this.sourceAssembly),
                this.backPropEnv.mongoTemplate, it),
             new RetryableBatchingCursor<>(
                     // No condition on remappedFrom attribute as above because merges
                     // can originate from a SS remapped in another assembly into SS remapped from the source assembly
                    where("_id").regex("^EVA3399_MERGED_${this.remappedAssembly}_.*"),
                    this.backPropEnv.mongoTemplate, it)]
        }
        opsOfImpactedSS.each{it.each{ it.each{List<SubmittedVariantOperationEntity> opsInBatch ->
            def impactedSSIds = opsInBatch.collect{it.accession}.toSet()
            impactedSSIds.addAll(opsInBatch.findAll{Objects.nonNull(it.mergedInto)}.collect{it.mergedInto}.toSet())
            [sveClass, dbsnpSveClass].each{ssClass ->
                this.backPropEnv.mongoTemplate.updateMulti(query(where("seq").is(this.sourceAssembly)
                        .and("accession").in(impactedSSIds)), new Update().unset("backPropRS"), ssClass)
            }
            def changedSS = [sveClass, dbsnpSveClass].collect{ssClass ->
                this.backPropEnv.mongoTemplate.find(query(where("seq").is(this.sourceAssembly)
                        .and("accession").in(impactedSSIds)), ssClass)
            }.flatten()
            backPropagatedRSWriter.write(changedSS)
        }}}
    }
}
