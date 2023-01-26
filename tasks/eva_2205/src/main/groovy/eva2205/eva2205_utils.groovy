package eva2205

import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment

import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

// Given some operations on SS, get the entire merge history for all the SS involved in those operations
static List getSSHistoryInvolvedInRSMerges(EVADatabaseEnvironment dbEnvToUse, String assembly, List svoeOps) {
    def hashesToLookFor = svoeOps.collect{it.inactiveObjects[0].hashedMessage}
    def ssMergeHistory = [svoeClass, dbsnpSvoeClass].collect{dbEnvToUse.mongoTemplate.find(query(where("inactiveObjects.seq").is(
            assembly).and("inactiveObjects.hashedMessage").in(hashesToLookFor).and("eventType").is(EventType.UPDATED).and(
            "reason").regex(".* was merged into .*")), it)}.flatten()
    return ssMergeHistory
}

// Get a chronological merge chain (oldest to most recent) given a set of RS merge operations that SS went through.
static List getChronologicalMergeChain(List svoeOpsForAGivenSS) {
    def opWithOriginalSS = null
    def getMatchedComponent = {op, component ->
        def matcher = op.reason =~ /Original rs(?<source>\d+) was merged into rs(?<destination>\d+)./
        if(matcher.matches()) {
            return matcher.group(component)
        }
        throw new IllegalArgumentException("Could not locate ${component} RS for ${op}!")
    }
    def opsGroupedByRSSource = svoeOpsForAGivenSS.groupBy {getMatchedComponent(it, "source")}
    def opsGroupedByRSDest = svoeOpsForAGivenSS.groupBy {getMatchedComponent(it, "destination")}
    opsGroupedByRSSource.each {srcRS, opsWithSrcRS ->
        if(!opsGroupedByRSDest.containsKey(srcRS)) {
            opWithOriginalSS = opsWithSrcRS[0]
        }
    }
    if (Objects.isNull(opWithOriginalSS)) throw new IllegalArgumentException("Could not locate source RS for the SS " +
            "with the following operations: \n" + svoeOpsForAGivenSS.join("\n"))
    def mergeChain = Arrays.asList(opWithOriginalSS)
    def destRSIDToLookFor = getMatchedComponent(opWithOriginalSS, "destination")
    while(true) {
        if(!opsGroupedByRSSource.containsKey(destRSIDToLookFor)) break
        def nextOpInChain = opsGroupedByRSSource.get(destRSIDToLookFor)[0]
        mergeChain += nextOpInChain
        destRSIDToLookFor = getMatchedComponent(nextOpInChain, "destination")
    }
    return mergeChain
}

// Given a chronological merge chain, get the most recent record without mapWeight
static def mostRecentNonMapWtSSRecord(List mergeChain) {
    return mergeChain.reverse().findAll{Objects.isNull(it.inactiveObjects[0].mapWeight)}[0]
}

// Given a chronological merge chain, get the most recent record with mapWeight
static def mostRecentMapWtSSRecord(List mergeChain) {
    return mergeChain.reverse().findAll{Objects.nonNull(it.inactiveObjects[0].mapWeight)}[0]
}
