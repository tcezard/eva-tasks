package src

@Grab(group='org.codehaus.groovy', module='groovy-cli-commons', version='3.0.10')

import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Criteria
import src.EVADatabaseEnvironment
import src.EVADataSet
import src.DuplicateSSCategory
import uk.ac.ebi.eva.accession.core.batch.io.SubmittedVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.deprecate.Application

import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where

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
def evaSSAccessionsLowerBound = destinationEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class)
def evaRSAccessionsLowerBound = destinationEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class)

def logger = LoggerFactory.getLogger(this.class)

List<String> assembliesToDeprecate =
        [SubmittedVariantEntity.class, DbsnpSubmittedVariantEntity.class].collect{collection ->
            sourceEnv.mongoTemplate.findDistinct(query(new Criteria()), "seq", collection, String.class)
        }.flatten().unique() as List<String>

assembliesToDeprecate.each {assembly ->
    def ssDeprecationWriter =
            new SubmittedVariantDeprecationWriter(assembly, destinationEnv.mongoTemplate,
                    destinationEnv.submittedVariantAccessioningService,
                    destinationEnv.clusteredVariantAccessioningService, evaSSAccessionsLowerBound, evaRSAccessionsLowerBound,
                    "EVA2936",
                    "Remapped variant had the same SS ID as the one imported from dbSNP. See EVA-2936.")
    logger.info("Deprecating impacted SS in assembly ${assembly}...")
    def deprecableSSDataSet = new EVADataSet(where("seq").is(assembly), sourceEnv.mongoTemplate,
            DuplicateSSCategory.class)
    deprecableSSDataSet.each{duplicateSSEntries ->
        List<Long> ssIDs = duplicateSSEntries.collect{it.accession}
        def svesToDeprecate = sourceEnv.submittedVariantAccessioningService.getAllActiveByAssemblyAndAccessionIn(assembly, ssIDs)
                .findAll{Objects.nonNull(it.getData().getRemappedFrom())}
                .collect{new SubmittedVariantEntity(it.getAccession(), it.getHash(), it.getData(), it.getVersion())}
        ssDeprecationWriter.write(svesToDeprecate)
    }
}
