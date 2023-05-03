package eva3238

import eva3238.DeprecateRemappedAllelesMatchFalseSS
import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.accession.deprecate.Application

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script deprecates variants remapped from variants having allelesMatch false in dbsnpSVE collection
// but had duplicate SS IDs in the remapped assembly
// So, remapping was carried out for the impacted variants in source assembly
// and the resultant remapped variants were placed in the DEV environment
// Since these resultant variants unambiguously show the remapped hashes, they can be readily deprecated from PROD
def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file to use for deprecation", required: true)
cli.devPropertiesFile(args: 1, "Development properties file to use", required: true)
cli.remappedAssembly(args: 1, "Remapped assembly", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class,
        ["parameters.assemblyAccession": options.remappedAssembly])
def devEnv = createFromSpringContext(options.devPropertiesFile, Application.class,
        ["parameters.assemblyAccession": options.remappedAssembly])
def ssHashesToDeprecate = devEnv.mongoTemplate.find(query(where("remappedFrom").exists(true)),
        dbsnpSveClass).collect{it.hashedMessage}
def svesToDeprecateInProd = [sveClass, dbsnpSveClass].collect {collectionClass ->
    prodEnv.mongoTemplate.find(query(where("_id").in(ssHashesToDeprecate)), collectionClass)
}.flatten()
DeprecateRemappedAllelesMatchFalseSS.deprecateSS(prodEnv, svesToDeprecateInProd)
