package eva3372

import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.accession.core.GenericApplication
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*


// This script collects pending splits and merges from the given assembly in PROD.
// If DEV environment is provided, will store them in
// eva3372_accession_sharded.rsIdsWithCveHash and cveHashesWithRsId collections.
// Otherwise will output intermediate results to file.
def cli = new CliBuilder()
cli.assemblyAccession(args: 1, "Assembly to analyze", required: true)
cli.prodPropertiesFile(args: 1, "Production properties file to use for database connection", required: true)
cli.devPropertiesFile(args: 1, "Development properties file to use for database connection", required: false)
cli.idAndHashFile(args: 1, "File to output IDs and hashes", required: false)

def options = cli.parse(args)
if (!options) {
    cli.usage()
    throw new Exception("Invalid command line options provided!")
}


// this is equivalent to if __name__ == '__main__' in Python
if (this.getClass().getName().equals('eva3372.eva3372_collect_rs_splits_and_merges')) {
    def prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
    if (options.devPropertiesFile) {
        def devEnv = createFromSpringContext(options.devPropertiesFile, GenericApplication.class)
        new CollectSplitsAndMerges(options.assemblyAccession, prodEnv, devEnv).collectAndCount()
    } else if (options.idAndHashFile) {
        new CollectIdsAndHashesToFile(options.assemblyAccession, prodEnv, options.idAndHashFile).exportRsIdsAndHashes()
    } else {
        throw new IllegalArgumentException("Must have either dev properties file or ID & hash file")
    }

}
