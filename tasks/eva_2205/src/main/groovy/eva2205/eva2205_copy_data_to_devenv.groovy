package eva2205

import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import groovy.cli.picocli.CliBuilder

def cli = new CliBuilder()
cli.propertiesFile(args:1, "Properties file to use for remediation", required: true)
cli.assembly(args:1, "Assembly to remediate", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
def devEnv = createFromSpringContext(options.devPropertiesFile, GenericApplication.class)
def assemblyAttributeInCollections = ["clusteredVariantEntity": "asm",
                                      "clusteredVariantOperationEntity": "inactiveObjects.asm", "dbsnpClusteredVariantEntity": "asm",
                                      "dbsnpClusteredVariantOperationEntity": "inactiveObjects.asm",
                                      "dbsnpSubmittedVariantEntity": "seq",
                                      "dbsnpSubmittedVariantOperationEntity": "inactiveObjects.seq",
                                      "submittedVariantEntity": "seq",
                                      "submittedVariantOperationEntity": "inactiveObjects.seq"]
def collectionClassForCollections = ["clusteredVariantEntity": cveClass,
                                     "clusteredVariantOperationEntity": cvoeClass, "dbsnpClusteredVariantEntity": dbsnpCveClass,
                                     "dbsnpClusteredVariantOperationEntity": dbsnpCvoeClass,
                                     "dbsnpSubmittedVariantEntity": dbsnpSveClass,
                                     "dbsnpSubmittedVariantOperationEntity": dbsnpSvoeClass,
                                     "submittedVariantEntity": sveClass,
                                     "submittedVariantOperationEntity": svoeClass]
def assembliesToCopy = ["GCA_000001515.4", "GCA_000001635.4"]
assembliesToCopy.each {assembly -> assemblyAttributeInCollections.each {collectionName, assemblyAttribute ->
    def collectionClass = collectionClassForCollections[collectionName]
    def cursorFromSourceCollection = new RetryableBatchingCursor(where(assemblyAttributeInCollections[collectionName]).
            is(assembly), prodEnv.mongoTemplate, collectionClass)
    cursorFromSourceCollection.each { entriesFromSource ->
        devEnv.bulkInsertIgnoreDuplicates(entriesFromSource, collectionClass)
    }
}}
