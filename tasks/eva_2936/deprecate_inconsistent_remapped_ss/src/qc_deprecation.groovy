package src

@Grab(group='org.codehaus.groovy', module='groovy-cli-commons', version='3.0.10')

import src.EVADatabaseEnvironment
import src.EVADataSet
import src.DuplicateSSCategory
import org.slf4j.LoggerFactory
import uk.ac.ebi.eva.accession.deprecate.Application

def cli = new CliBuilder(usage: 'deprecate_eva2936_variants <source properties file> <destination properties file>')
def options = cli.parse(args)
if (options.arguments().size() != 2) {
    cli.usage()
    return
}

// Both can be specified as the same file
// in which case variants to be deprecated will be read from/output written in the same environment
def sourcePropertiesFile = options.arguments()[0]
def destinationPropertiesFile = options.arguments()[1]
def sourceEnv = EVADatabaseEnvironment.createFromSpringContext(sourcePropertiesFile, Application.class)
def destinationEnv = (destinationPropertiesFile == sourcePropertiesFile)? sourceEnv: EVADatabaseEnvironment.createFromSpringContext(destinationPropertiesFile, Application.class)

def logger = LoggerFactory.getLogger(this.class)
def deprecableSSDataset = new EVADataSet(null, sourceEnv.mongoTemplate, DuplicateSSCategory.class)

// Ensure that there are no remaining accessions
// in the three categories mentioned here: https://docs.google.com/spreadsheets/d/1LZbUtamFtQH12SVTfCfYOjcIHWL1OOC_kYnieRLR9UA/edit#rangeid=1722186521
deprecableSSDataset.each{it.groupBy {it.getSeq()}.each{String assembly, List<DuplicateSSCategory> dupSSEntries ->
    List<Long> ssIDs = dupSSEntries.collect{it.getAccession()}
    destinationEnv.submittedVariantAccessioningService.getAllActiveByAssemblyAndAccessionIn(assembly, ssIDs)
            .groupBy {it.getAccession()}.findAll{accession, sves ->
        // There should not be any accessions from multiple source assemblies
        return sves.collect{it.getData().getRemappedFrom()}.unique().size() > 1
    }.each{accession, svesStillInDeprecableCategories ->
        logger.error("Remapped variants with ss${accession} in ${assembly} should have been deprecated!!")}
}}
