package eva3399

import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.eva.accession.clustering.batch.io.BackPropagatedRSWriter
import uk.ac.ebi.eva.accession.clustering.configuration.BeanNames
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
        def impactedSS = [sveClass, dbsnpSveClass].collect{new RetryableBatchingCursor<>(
                where("seq").is(this.remappedAssembly).and("remappedFrom").is(this.sourceAssembly)
                        .and("rs").exists(true), this.backPropEnv.mongoTemplate, it)}
        impactedSS.each{it.each{ssInBatch ->
            def impactedSSIds = ssInBatch.collect{it.accession}.toSet()
            [sveClass, dbsnpSveClass].each{ssClass ->
                this.backPropEnv.mongoTemplate.updateMulti(query(where("seq").is(this.sourceAssembly)
                        .and("accession").in(impactedSSIds)), new Update().unset("backPropRS"), ssClass)
            }
            def changedSS = [sveClass, dbsnpSveClass].collect{ssClass ->
                this.backPropEnv.mongoTemplate.find(query(where("seq").is(this.sourceAssembly)
                        .and("accession").in(impactedSSIds)), ssClass)
            }.flatten()
            backPropagatedRSWriter.write(changedSS)
        }}
    }
}
