package eva3450

import groovy.cli.picocli.CliBuilder
import org.springframework.boot.SpringApplication
import org.springframework.boot.autoconfigure.SpringBootApplication
import org.springframework.context.ConfigurableApplicationContext
import uk.ac.ebi.eva.accession.core.GenericApplication
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*


// This script collects pending splits and merges from the given assembly in PROD.
// If DEV environment is provided, will store them in
// eva3372_accession_sharded.rsIdsWithCveHash and cveHashesWithRsId collections.
// Otherwise will output intermediate results to file.
def cli = new CliBuilder()
cli.rsMergeCandidateFile(args: 1, "file containing the RS merge candidates (Single Hash multiple rs)", required: true)
cli.assemblyAccession(args: 1, "Assembly accession where the submitted variant should be found", required: true)
cli.devPropertiesFile(args: 1, "Properties file to use for the dev database connection", required: true)
cli.prodPropertiesFile(args: 1, "Properties file to use for the prod database connection", required: true)
cli.batchSize(args: 1, "Properties file to use for database connection", required: false, defaultValue: "1000")


def options = cli.parse(args)
if (!options) {
    cli.usage()
    throw new Exception("Invalid command line options provided!")
}

// this is equivalent to if __name__ == '__main__' in Python
if (this.getClass().getName().equals('eva3450.eva3450_create_merge_candidates')) {
    def devEnv = createFromSpringContext(options.devPropertiesFile, GenericApplication.class)
    def prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
    new CreateMergeCandidates(options.assemblyAccession, options.rsMergeCandidateFile, prodEnv, devEnv, options.batchSize.toInteger()).process()
}

