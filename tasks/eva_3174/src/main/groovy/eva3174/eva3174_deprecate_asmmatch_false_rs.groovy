package eva3174

import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory

import uk.ac.ebi.eva.accession.clustering.metric.ClusteringMetric
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.batch.io.ClusteredVariantDeprecationWriter
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.EVAObjectModelUtils
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.metrics.metric.MetricCompute

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

def cli = new CliBuilder()
cli.propertiesFile(args:1, "Properties file to use for deprecation", required: true)
cli.assemblyToDeprecate(args:1, "Assembly where RS corresponding to asmMatch=false should be deprecated", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

class DeprecateRSWithAsmMatchFalseLocus {
    EVADatabaseEnvironment dbEnv
    MetricCompute metricCompute
    String assembly
    Long dbsnpRSIDUpperBound
    static def logger = LoggerFactory.getLogger(GenericApplication.class)

    DeprecateRSWithAsmMatchFalseLocus(EVADatabaseEnvironment dbEnv, String assembly) {
        this.dbEnv = dbEnv
        this.assembly = assembly
        this.dbsnpRSIDUpperBound = dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class)
        this.metricCompute = dbEnv.springApplicationContext.getBean(MetricCompute.class)
    }

    def deprecateRS() {
        def asmMatchCursorDbsnpSve = [dbsnpSvoeClass, svoeClass].collect{new RetryableBatchingCursor<>(
                where("inactiveObjects.seq").is(this.assembly).and("_id").regex("SS_DEPRECATED_EVA3174_.*"),
                this.dbEnv.mongoTemplate, it)}.flatten()
        asmMatchCursorDbsnpSve.each {it.each{deprecatedAsmMatchOps ->
            def cveHashesToDeprecate = deprecatedAsmMatchOps.collect{
                EVAObjectModelUtils.getClusteredVariantHash(it.inactiveObjects[0])}
            def cvesToDeprecate = [cveClass, dbsnpCveClass].collect{
                dbEnv.mongoTemplate.find(query(where("_id").in(cveHashesToDeprecate)), it)}.flatten()
            logger.info("Deprecating ${cvesToDeprecate.size()} CVEs in assembly ${this.assembly}...")
            def cvDeprecationWriter = new ClusteredVariantDeprecationWriter(this.assembly,
                    dbEnv.mongoTemplate, dbEnv.submittedVariantAccessioningService, dbsnpRSIDUpperBound,
                    "EVA3174", "Variant deprecated due to allelesMatch=false")
            // By the design of the write, this will only succeed if these RS IDs are orphaned i.e., there are no SS that have these RS IDs
            cvDeprecationWriter.write(cvesToDeprecate)
            this.metricCompute.addCount(ClusteringMetric.CLUSTERED_VARIANTS_DEPRECATED,
                    cvDeprecationWriter.numDeprecatedEntities)
            this.metricCompute.saveMetricsCountsInDB()
        }}
    }
}

EVADatabaseEnvironment dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assemblyToDeprecate])
def deprecateObj = new DeprecateRSWithAsmMatchFalseLocus(dbEnv, options.assemblyToDeprecate)
deprecateObj.deprecateRS()
