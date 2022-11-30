package uk.ac.ebi.eva.eva3100

import groovy.cli.picocli.CliBuilder
import org.springframework.boot.SpringApplication
import org.springframework.data.mongodb.core.query.Update
import org.slf4j.LoggerFactory
import uk.ac.ebi.eva.accession.core.GenericApplication

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where

// This script restores the correct allelesMatch attribute inside dbsnpSVOE operations
// with SS_SPLIT_FROM involving SS with questionable allelesMatch

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning", required: true)
cli.devPropertiesFile(args:1, "Development properties file for accessioning", required: true)
cli.allelesMatchSplitSSLog(args:1, "ss_split_from_ss_with_alleles_match_false file from EVA-3100", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
// Environment pointing to pre-EVA2861 split
def devEnv = createFromSpringContext(options.devPropertiesFile, GenericApplication.class)
def scriptLogger = LoggerFactory.getLogger(SpringApplication.class)

// dbsnp SVOE IDs which report split from dbSNP SS with allelesMatch in EVA-3100
def impactedIDs = ["bash", "-c", "cut -d, -f1 ${options.allelesMatchSplitSSLog}"].execute().text.split("\n")
//436,667
println(impactedIDs.size())

// Check if any of the SS  were split from
def numIDsScanned = 0
impactedIDs.collate(1000).each {dbsnpSvoeIdsToLookFor ->
    def impactedSvoes = prodEnv.mongoTemplate.find(query(where("_id").in(dbsnpSvoeIdsToLookFor)), dbsnpSvoeClass)
    def hashesToFind = impactedSvoes.collect{it.inactiveObjects[0].hashedMessage}.toSet()
    numIDsScanned += hashesToFind.size()
    // get the SS record from dbsnpSVE before it was split
    def preSplitSS = devEnv.mongoTemplate.find(query(where("_id").in(hashesToFind)), dbsnpSveClass,
            "dbsnpSubmittedVariantEntity_before_eva2861_split")
    def matchingHashes = preSplitSS.collect { it.hashedMessage }.toSet()
    scriptLogger.info("${preSplitSS.size()} matching SS found...")
    (hashesToFind - matchingHashes).each { scriptLogger.warn("Could not find hash ${it}...") }
    // Did the original SS records before split have the allelesMatch flag?
    def matchingOriginalSSWithNoAllelesMatch = preSplitSS.findAll { !(it.isAllelesMatch()) }
    if (matchingOriginalSSWithNoAllelesMatch.size() > 0) {
        scriptLogger.warn(matchingOriginalSSWithNoAllelesMatch.size() + " SS have allelesMatch false!!")
        def hashesWithNoAllelesMatch = matchingOriginalSSWithNoAllelesMatch.collect { it.hashedMessage }.toSet()
        hashesWithNoAllelesMatch.each{scriptLogger.warn("SS hash ${it} has allelesMatch false!!")}
        def dbsnpSvoeIDsToUpdate = impactedSvoes.findAll { hashesWithNoAllelesMatch.contains(it.inactiveObjects[0].hashedMessage) }.collect { it.getId() }
        dbsnpSvoeIDsToUpdate.each{scriptLogger.warn("DbsnpSVOE with ID ${it} will be updated!!")}
        prodEnv.mongoTemplate.updateMulti(query(where("_id").in(dbsnpSvoeIDsToUpdate)),
                Update.update("inactiveObjects[0].allelesMatch", false), dbsnpSvoeClass)
    }
    scriptLogger.info("${numIDsScanned} dbsnpSVE hashes scanned so far...")
}

