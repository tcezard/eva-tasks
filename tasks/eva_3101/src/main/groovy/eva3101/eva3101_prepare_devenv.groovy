package eva3101

import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

// This script populates the development environment test database with
// variants having allelesMatch attribute false in dbsnpSVE collection
def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file for accessioning", required: true)
cli.devPropertiesFile(args: 1, "Development properties file for accessioning", required: true)
cli.options.assemblyToDeprecate(args: 1, "Assembly to be deprecated", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class)
def devEnv = createFromSpringContext(options.devPropertiesFile, Application.class)

// Transfer data to DEV for an assembly
def ssEntryBatches = new RetryableBatchingCursor(where("seq").is(options.assemblyToDeprecate).and("allelesMatch").exists(true),
        prodEnv.mongoTemplate, dbsnpSveClass)
ssEntryBatches.each { ssEntries ->
    devEnv.bulkInsertIgnoreDuplicates(ssEntries, dbsnpSveClass)
    List<Long> correspondingRSIDs = ssEntries.collect { it.clusteredVariantAccession }.findAll {Objects.nonNull(it)}
    List<String> correspondingSSHashes = ssEntries.collect { it.hashedMessage }
    def correspondingRSEntries = new RetryableBatchingCursor(where("asm").is(options.assemblyToDeprecate)
            .and("accession").in(correspondingRSIDs), prodEnv.mongoTemplate,
            dbsnpCveClass).collect()
    def otherSSEntriesWithSameRS = new RetryableBatchingCursor(where("seq").is(options.assemblyToDeprecate)
            .and("rs").in(correspondingRSIDs).and("_id").nin(correspondingSSHashes),
            prodEnv.mongoTemplate, dbsnpSveClass).collect()
    devEnv.bulkInsertIgnoreDuplicates(correspondingRSEntries, dbsnpCveClass)
    devEnv.bulkInsertIgnoreDuplicates(otherSSEntriesWithSameRS, dbsnpSveClass)
}
