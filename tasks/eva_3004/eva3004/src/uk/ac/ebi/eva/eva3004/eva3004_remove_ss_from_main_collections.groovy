// This script removes the SS shelved in EVA-3004 from the main dbsnpSVE and SVE collections

package uk.ac.ebi.eva.eva3004

import groovy.cli.picocli.CliBuilder
import org.springframework.data.mongodb.core.query.Criteria
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.eva3004.EVACursor
import uk.ac.ebi.eva.eva3004.EVALoggingUtils

import static uk.ac.ebi.eva.eva3004.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning",  required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class)
def scriptLogger = EVALoggingUtils.getLogger(this.class)

def dbsnpSveShelvedCollections = ["eva2979_dbsnpSubmittedVariantEntity", "eva2706_shelved_dbsnpsve_multi_position_rs",
                                  "eva2915_shelved_dbsnpsve_dups_wo_mapwt", "eva2950_dbsnpSubmittedVariantEntity",
                                  "eva2959_shelved_dbsnpsve"]
def sveShelvedCollections = ["eva2706_shelved_sve_multi_position_rs", "eva2915_shelved_sve_split_from_mapwt_ss",
                             "eva2950_submittedVariantEntity", "eva2959_shelved_sve"]

dbsnpSveShelvedCollections.each{collectionName ->
    def dbsnpSveCursorFromShelvedCollection = new EVACursor(new Criteria(), prodEnv.mongoTemplate, dbsnpSveClass,
            pageSize = 1000, collectionName)
    def numRemovedDocuments = 0
    dbsnpSveCursorFromShelvedCollection.each{dbsnpSvesInBatch ->
        def removedDocs = prodEnv.mongoTemplate.findAllAndRemove(
                query(where("_id").in(dbsnpSvesInBatch.collect{it.getId()})), dbsnpSveClass)
        numRemovedDocuments += removedDocs.size()
        scriptLogger.info("Removed ${numRemovedDocuments} so far from dbsnpSubmittedVariantEntity using " +
                "${collectionName} as reference!!")
    }
}

sveShelvedCollections.each{collectionName ->
    def sveCursorFromShelvedCollection = new EVACursor(new Criteria(), prodEnv.mongoTemplate, sveClass,
            pageSize = 1000, collectionName)
    def numRemovedDocuments = 0
    sveCursorFromShelvedCollection.each{svesInBatch ->
        def removedDocs = prodEnv.mongoTemplate.findAllAndRemove(
                query(where("_id").in(svesInBatch.collect{it.getId()})), sveClass)
        numRemovedDocuments += removedDocs.size()
        scriptLogger.info("Removed ${numRemovedDocuments} so far from submittedVariantEntity using " +
                "${collectionName} as reference!!")
    }
}
