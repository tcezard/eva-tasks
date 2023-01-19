package eva2205

import org.bson.types.ObjectId
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment

import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

// Given some operations on SS, get the entire merge history for all the SS involved in those operations
static List getSSHistoryInvolvedInRSMerges(EVADatabaseEnvironment dbEnvToUse, String assembly, List svoeOps) {
    def accessionsToLookFor = svoeOps.collect{it.accession}
    def validAccessionHashCombos = svoeOps.collect { it.accession + "_" +
            it.inactiveObjects[0].hashedMessage }.toSet()
    def ssMergeHistory = [svoeClass, dbsnpSvoeClass].collect{dbEnvToUse.mongoTemplate.find(query(where("inactiveObjects.seq").is(
            assembly).and("accession").in(accessionsToLookFor).and("eventType").is(EventType.UPDATED).and(
            "reason").regex(".* was merged into .*")), it)}.flatten().findAll{op ->
        validAccessionHashCombos.contains(op.accession + "_" + op.inactiveObjects[0].hashedMessage)}
    return ssMergeHistory
}

// Get original SS record given a set of operations on that SS. Do this in-memory because doing this for thousands of SS in-database could get expensive!
static List getOriginalSSRecords(List svoeOpsForAGivenSS) {
    def originalSSRecords = null
    def getMatchedComponent = {op, component ->
        def matcher = op.reason =~ /Original rs(?<source>\d+) was merged into rs(?<destination>\d+)./
        if(matcher.matches()) {
            return matcher.group(component)
        }
        throw new IllegalArgumentException("Could not locate ${component} RS for ${op}!")
    }
    def opsGroupedByRSSource = svoeOpsForAGivenSS.groupBy {
        getMatchedComponent(it, "source")
    }
    def opsGroupedByRSDest = svoeOpsForAGivenSS.groupBy {
        getMatchedComponent(it, "destination")
    }
    opsGroupedByRSSource.each {srcRS, opsWithSrcRS ->
        if(!opsGroupedByRSDest.containsKey(srcRS)) {
            originalSSRecords = opsWithSrcRS.collect{
                it.inactiveObjects[0].toSubmittedVariantEntity()
            }
        }
    }
    if (Objects.isNull(originalSSRecords)) throw new IllegalArgumentException("Could not locate source RS for the SS " +
            "with the following operations: \n" + svoeOpsForAGivenSS.join("\n"))
    return originalSSRecords
}

def svoeOps = dbEnv.mongoTemplate.find(query(where("_id").is(new ObjectId("5d3a4acc7452826c5d3fbbd9"))), dbsnpSvoeClass)
println(getOriginalSSRecords(getSSHistoryInvolvedInRSMerges(dbEnv, options.assembly, svoeOps)))