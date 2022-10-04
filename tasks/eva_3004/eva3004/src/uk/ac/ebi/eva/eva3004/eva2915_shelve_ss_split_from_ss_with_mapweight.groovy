/***
 * This script shelves SS entries in EVA which were split from dbsnp SS which have mapping weights
 */

package uk.ac.ebi.eva.eva3004

import groovy.cli.picocli.CliBuilder
import org.springframework.data.mongodb.core.query.Query
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.eva3004.EVACursor

import static uk.ac.ebi.eva.eva3004.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, longOpt: "prod-props-file", "Production properties file for accessioning",  required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}
def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class)

def allSveSplitOps = new EVACursor(where("_id").regex("SS_SPLIT_.*"), prodEnv.mongoTemplate, dbsnpSvoeClass)
def shelvedCollectionSve = "eva2915_shelved_sve_split_from_mapwt_ss"
def totalSvesShelved = 0
allSveSplitOps.each{ops -> ops.groupBy{it.inactiveObjects[0].referenceSequenceAccession}.each{assembly, sveSplitOpsInBatch ->
    def splitOpsGroupedByAccession = sveSplitOpsInBatch.groupBy{it.accession}
    def sveAccessionsToLookup = splitOpsGroupedByAccession.keySet()
    def dbsnpSvesWithMapWeight = prodEnv.mongoTemplate.find(query(where("seq").is(assembly).and("accession").in(sveAccessionsToLookup)
            .and("mapWeight").exists(true)), dbsnpSveClass).findAll{it.mapWeight > 1}.collect{it.accession}
    def evaSveHashesSplitFromDbsnpSvesWithMapWeight= dbsnpSvesWithMapWeight.collect{splitOpsGroupedByAccession.get(it).collect{it.inactiveObjects[0].hashedMessage}}.flatten()
    def evaSvesToShelve =
            prodEnv.mongoTemplate.find(query(where("_id").in(evaSveHashesSplitFromDbsnpSvesWithMapWeight)), sveClass)
    println("Shelving ${evaSvesToShelve.size()} EVA SS...")
    prodEnv.bulkInsertIgnoreDuplicates(evaSvesToShelve, sveClass, shelvedCollectionSve)
    totalSvesShelved += evaSvesToShelve.size()
}}
// 280,606
println(prodEnv.mongoTemplate.count(new Query(), sveClass, shelvedCollectionSve))
