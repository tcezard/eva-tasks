package eva2205

import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory

import uk.ac.ebi.eva.accession.clustering.metric.ClusteringMetric
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.batch.io.ClusteredVariantDeprecationWriter
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.metrics.metric.MetricCompute

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

def cli = new CliBuilder()
cli.propertiesFile(args:1, "Properties file to use for remediation", required: true)
cli.assembly(args:1, "Assembly to remediate", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

class DeprecateOrphanedRS {
    EVADatabaseEnvironment dbEnv
    MetricCompute metricCompute
    String assembly
    Long dbsnpRSIDUpperBound
    static def logger = LoggerFactory.getLogger(GenericApplication.class)

    DeprecateOrphanedRS(EVADatabaseEnvironment dbEnv, String assembly) {
        this.dbEnv = dbEnv
        this.assembly = assembly
        this.dbsnpRSIDUpperBound = dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class)
        this.metricCompute = dbEnv.springApplicationContext.getBean(MetricCompute.class)
    }

    def deprecateOrphanedCves() {
        def mapWtCursorDbsnpCve =
                new RetryableBatchingCursor<>(where("asm").is(this.assembly).and("mapWeight").exists(true),
                        this.dbEnv.mongoTemplate, dbsnpCveClass)
        mapWtCursorDbsnpCve.each {cvesToDeprecate ->
            logger.info("Deprecating ${cvesToDeprecate.size()} CVEs in assembly ${this.assembly}...")
            def cvDeprecationWriter = new ClusteredVariantDeprecationWriter(this.assembly,
                    dbEnv.mongoTemplate, dbEnv.submittedVariantAccessioningService, dbsnpRSIDUpperBound,
                    "EVA2205", "Variant deprecated due to mapWeight > 1")
            // By the design of the write, this will only succeed if these RS IDs are orphaned i.e., there are no SS that have these RS IDs
            cvDeprecationWriter.write(cvesToDeprecate)
            this.metricCompute.addCount(ClusteringMetric.CLUSTERED_VARIANTS_DEPRECATED,
                    cvDeprecationWriter.numDeprecatedEntities)
            this.metricCompute.saveMetricsCountsInDB()
        }
    }
}

EVADatabaseEnvironment dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assembly])
def deprecateOrphanedRSObj = new DeprecateOrphanedRS(dbEnv, options.assembly)
deprecateOrphanedRSObj.deprecateOrphanedCves()
