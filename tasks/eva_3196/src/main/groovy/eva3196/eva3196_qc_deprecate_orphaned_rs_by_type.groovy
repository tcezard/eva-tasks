package eva3196

import org.apache.commons.lang3.tuple.ImmutablePair
import org.springframework.data.mongodb.core.query.Criteria
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// QC - EVA3196
// Ensure no SS is associated with deprecated RS locus
// No error messages were printed below
def cvOpCursor = [cvoeClass, dbsnpCvoeClass].collect{opClass ->
    new RetryableBatchingCursor<>(where("_id").regex("RS_DEPRECATED_EVA3196_.*"), dbEnv.mongoTemplate, opClass)
}
int numEntriesScanned = 0
cvOpCursor.each {it.each{deprecationOps ->
    numEntriesScanned += deprecationOps.size()
    println("Scanned ${numEntriesScanned} ops entries so far...")
    def deprecatedRSHashesAndIDs = deprecationOps.collect{
        new ImmutablePair<>(it.inactiveObjects[0].hashedMessage, it.accession)}.toSet()
    def svesWithDeprecatedRS = [sveClass, dbsnpSveClass].collect{collectionClass ->
        dbEnv.mongoTemplate.find(query(where("seq").is(options.assemblyToDeprecate)
                .and("rs").in(deprecatedRSHashesAndIDs.collect{it.right})),
                collectionClass)
                .findAll{deprecatedRSHashesAndIDs.contains(
                        new ImmutablePair<>(EVAObjectModelUtils.getClusteredVariantHash(it) ,
                                it.clusteredVariantAccession))}
    }.flatten()
    svesWithDeprecatedRS.each{println("ERROR: This submitted variant still has RS locus " +
            "of a deprecated RS: ${it}")}
}}

// Ensure that SS are still associated with any RS left behind undeprecated
// No error messages were printed below
def cveCursor = [cveClass, dbsnpCveClass].collect{
    new RetryableBatchingCursor<>(new Criteria(), dbEnv.mongoTemplate, it)
}
numEntriesScanned = 0
cveCursor.each {it.each{cves ->
    numEntriesScanned += cves.size()
    println("Scanned ${numEntriesScanned} CVE entries so far...")
    def undeprecatedRSHashesAndIDsFromCVE = cves.collect{
        new ImmutablePair<>(it.hashedMessage, it.accession)}.toSet()
    def undeprecatedRSHashesAndIDsFromSVE = [sveClass, dbsnpSveClass].collect{collectionClass ->
        dbEnv.mongoTemplate.find(query(where("seq").is(options.assemblyToDeprecate)
                .and("rs").in(undeprecatedRSHashesAndIDsFromCVE.collect{it.right})),
                collectionClass)
                .findAll{undeprecatedRSHashesAndIDsFromCVE.contains(
                        new ImmutablePair<>(EVAObjectModelUtils.getClusteredVariantHash(it),
                                it.clusteredVariantAccession))}
    }.flatten().collect{new ImmutablePair<>(EVAObjectModelUtils.getClusteredVariantHash(it) ,
            it.clusteredVariantAccession)}.toSet()
    (undeprecatedRSHashesAndIDsFromCVE - undeprecatedRSHashesAndIDsFromSVE).each{
        println("ERROR: This RS locus and accession of an undeprecated RS has no associated SS: ${it}")}
}}
