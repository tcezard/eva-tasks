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
    /** There are scenarios (see example below) where there could be circular merges leaving no way to find a definitive original SS record!!
     * In those cases, use Submitted Variant createdDate (NOT the createdDate for the merge operation which was unfortunately not available for old records)
     * as a best-effort proxy (possibly not accurate!) to identify the original SS record
     * DbsnpSubmittedVariantOperationEntity{id='5d554066d7d16e435e2b24ea', accession=316150381, reason='Original rs196349204 was merged into rs195124461.'
     DbsnpSubmittedVariantOperationEntity{id='5d554066d7d16e435e2b24ee', accession=315919449, reason='Original rs195124461 was merged into rs196349204.'
     */
    if (Objects.isNull(opWithOriginalSS)) {
        opWithOriginalSS = svoeOpsForAGivenSS.sort { op1, op2 ->
            op1.inactiveObjects[0].createdDate.compareTo(op2.inactiveObjects[0].createdDate)}[0]
    }
    def mergeChain = Arrays.asList(opWithOriginalSS)
    def destRSIDToLookFor = getMatchedComponent(opWithOriginalSS, "destination")
    def previouslyEncounteredRS = new HashSet<String>([getMatchedComponent(opWithOriginalSS, "source")])
    while(true) {
        if(!opsGroupedByRSSource.containsKey(destRSIDToLookFor)) break
        previouslyEncounteredRS.add(destRSIDToLookFor)
        def nextOpInChain = opsGroupedByRSSource.get(destRSIDToLookFor)[0]
        mergeChain += nextOpInChain
        destRSIDToLookFor = getMatchedComponent(nextOpInChain, "destination")
        // When encountering circular reference ex: rs3 -> rs2 -> rs1 -> rs3
        // break chain at rs1 whose destRSIDToLookFor will be rs3
        // See https://www.ebi.ac.uk/panda/jira/browse/EVA-3151?focusedId=418704&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-418704
        if(previouslyEncounteredRS.contains(destRSIDToLookFor)) break
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
