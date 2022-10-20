// This script restores the SS shelved in EVA-3004 into the main dbsnpSVE and SVE collections

package uk.ac.ebi.eva.eva3027

import groovy.cli.picocli.CliBuilder
import org.springframework.data.mongodb.core.query.Criteria
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.eva3027.EVACursor
import uk.ac.ebi.eva.eva3027.EVADatabaseEnvironment
import uk.ac.ebi.eva.eva3027.EVALoggingUtils

import static uk.ac.ebi.eva.eva3027.EVADatabaseEnvironment.*

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning", required: true)
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
    def numInsertedDocuments = 0
    dbsnpSveCursorFromShelvedCollection.each{dbsnpSvesInBatch ->
        def insertedCount = prodEnv.bulkInsertIgnoreDuplicates(dbsnpSvesInBatch, dbsnpSveClass)
        scriptLogger.info("Inserted ${numInsertedDocuments} so far from dbsnpSubmittedVariantEntity using " +
                "${collectionName} as reference!!")
        numInsertedDocuments += insertedCount
    }
}

sveShelvedCollections.each{collectionName ->
    def sveCursorFromShelvedCollection = new EVACursor(new Criteria(), prodEnv.mongoTemplate, sveClass,
            pageSize = 1000, collectionName)
    def numInsertedDocuments = 0
    sveCursorFromShelvedCollection.each{svesInBatch ->
        def insertedCount = prodEnv.bulkInsertIgnoreDuplicates(svesInBatch, sveClass)
        scriptLogger.info("Inserted ${numInsertedDocuments} so far from SubmittedVariantEntity using " +
                "${collectionName} as reference!!")
        numInsertedDocuments += insertedCount
    }
}
