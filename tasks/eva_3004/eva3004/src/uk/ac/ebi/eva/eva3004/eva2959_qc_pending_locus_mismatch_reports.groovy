// This is mostly a copy of qc_eva2959_prod. However, it only checks for pending RS/SS locus mismatches but does NOT shelve any entries.

package uk.ac.ebi.eva.eva3004

import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.eva3004.EVADatabaseEnvironment
import uk.ac.ebi.eva.eva3004.EVALoggingUtils

import static uk.ac.ebi.eva.eva3004.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, longOpt: "prod-props-file", "Production properties file for accessioning",  required: true)
cli.discordantRSDir(args:1, longOpt: "eva3004-discordant-rs-dir", "Directory containing logs of discordant RS " +
        "for each assembly ",  required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}
def (prodPropertiesFile, discordantRSDir) = [options.prodPropertiesFile, options.discordantRSDir]

def prodEnv = createFromSpringContext(prodPropertiesFile, Application.class)
def scriptLogger = EVALoggingUtils.getLogger(this.class)

// SS/RS locus mismatches reported by re-running scripts from EVA-2706 in
// Two categories:
// 1) More than one RS locus entry in the clustered variant for a given RS accession (not necessarily that these multiple loci are observed in the assigned SS)
// 2) More than one position reported for the same RS across multiple SS
def rsToCheckCategory1 = ["bash", "-c",
                    "ls -1 ${discordantRSDir}/*pending*errors.log"].execute().text.split("\n").collectEntries{
    [it.split("/").reverse()[0].split("_pending_errors")[0],
     ["bash", "-c", "grep \"cluster position\" ${it} | grep -o -E \"for rs[0-9]+\" " +
             "| grep -o -E \"[0-9]+\""].execute().text.trim().split("\n").findAll{!it.empty}.collect{it.toLong()}.toSet()
    ]}
def rsToCheckCategory2 = ["bash", "-c",
                    "ls -1 ${discordantRSDir}/*pending*errors.log"].execute().text.split("\n").collectEntries{
    [it.split("/").reverse()[0].split("_pending_errors")[0],
     ["bash", "-c", "grep \"positions found in submitted variants\" ${it} | grep -o -E \"for rs[0-9]+\" " +
             "| grep -o -E \"[0-9]+\""].execute().text.trim().split("\n").findAll{!it.empty}.collect{it.toLong()}.toSet()
    ]}
// 64,766
println(rsToCheckCategory1.values().flatten().size())
// 259
println(rsToCheckCategory2.values().flatten().size())

def shelvedCollectionDbsnpSve = "eva2706_shelved_dbsnpsve_multi_position_rs"
def shelvedCollectionSve = "eva2706_shelved_sve_multi_position_rs"
def shelveSSWithDiscordantRS = {String assembly, List<Long> rsIDs, String category ->
    List<SubmittedVariantEntity> allSvesWithRS = [sveClass, dbsnpSveClass].collect{prodEnv.mongoTemplate.find(query(where("seq").is(assembly)
            .and("rs").in(rsIDs)), it)}.flatten()
    // Ensure that the RS reported with multiple positions among the SS still holds true in the current SVE collections
    def allSvesToShelve = allSvesWithRS.groupBy{it.clusteredVariantAccession}.findAll{k, v ->
        v.collect{sve -> EVADatabaseEnvironment.toClusteredVariantEntity(sve).hashedMessage}.unique().size() > 1}.values().flatten()

    def dbsnpSvesToShelve = allSvesToShelve.findAll{it.accession < 5e9}
    def evaSvesToShelve = allSvesToShelve.findAll{it.accession >= 5e9}

    if (dbsnpSvesToShelve.size() > 0) {
        scriptLogger.error("Pending ${dbsnpSvesToShelve.size()} dbsnp SVEs in assembly ${assembly} in ${category} will be shelved!!")
        scriptLogger.error("${dbsnpSvesToShelve.collect{it.accession}}")
        prodEnv.bulkInsertIgnoreDuplicates(dbsnpSvesToShelve, dbsnpSveClass, shelvedCollectionDbsnpSve)
        prodEnv.mongoTemplate.findAllAndRemove(query(where("_id").in(dbsnpSvesToShelve.collect{it.id})), dbsnpSveClass)
    }
    if (evaSvesToShelve.size() > 0) {
        scriptLogger.error("Pending ${evaSvesToShelve.size()} EVA SVEs in assembly ${assembly} in ${category} will be shelved!!")
        scriptLogger.error("${evaSvesToShelve.collect{it.accession}}")
        prodEnv.bulkInsertIgnoreDuplicates(evaSvesToShelve, sveClass, shelvedCollectionSve)
        prodEnv.mongoTemplate.findAllAndRemove(query(where("_id").in(evaSvesToShelve.collect{it.id})), sveClass)
    }
    scriptLogger.info("Scanned ${rsIDs.size()} RS IDs for assembly ${assembly}...")
}

def category = "category 1"
rsToCheckCategory1.each{assembly, allRSIDs -> allRSIDs.collate(1000).each { rsIDs ->
    shelveSSWithDiscordantRS(assembly, rsIDs, category)
}}
category = "category 2"
rsToCheckCategory2.each{assembly, allRSIDs -> allRSIDs.collate(1000).each { rsIDs ->
    shelveSSWithDiscordantRS(assembly, rsIDs, category)
}}
