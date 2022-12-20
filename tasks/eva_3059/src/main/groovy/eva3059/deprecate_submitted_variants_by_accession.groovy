package eva3059

import groovy.cli.picocli.CliBuilder
import org.springframework.boot.SpringApplication
import org.slf4j.LoggerFactory
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.batch.io.SubmittedVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*


def cli = new CliBuilder()
cli.propertiesFile(args:1, "Properties file for accessioning", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def env = createFromSpringContext(options.propertiesFile, GenericApplication.class)
def logger = LoggerFactory.getLogger(SpringApplication.class)

def ssIdToDeprecate = 5015497294
def assemblyAccession = "GCA_000181335.3"

def sveWrapper = env.submittedVariantAccessioningService.getByAccession(ssIdToDeprecate)
def sveToDeprecate = new SubmittedVariantEntity(sveWrapper.getAccession(), sveWrapper.getHash(), sveWrapper.getData(), sveWrapper.getVersion())

logger.info("Deprecating ss${ssIdToDeprecate} in assembly ${assemblyAccession}...")

def writer = new SubmittedVariantDeprecationWriter(
        assemblyAccession,
        env.mongoTemplate,
        env.submittedVariantAccessioningService,
        env.clusteredVariantAccessioningService,
        env.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class),
        env.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
        "EVA3059",
        "Invalid variant deprecated during EVA3059"
)
writer.write([sveToDeprecate])
