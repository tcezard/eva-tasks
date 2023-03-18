package eva3196

import groovy.cli.picocli.CliBuilder
import org.apache.commons.lang3.tuple.ImmutablePair
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.clustering.metric.ClusteringMetric
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.accession.core.batch.io.ClusteredVariantDeprecationWriter
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.metrics.metric.MetricCompute

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script deprecates orphaned RS by type
// i.e., RS with a given RS locus not associated to any SS with that RS locus
def cli = new CliBuilder()
cli.propertiesFile(args: 1, "Production properties file to use for deprecation", required: true)
cli.assemblyToDeprecate(args: 1, "Assembly where deprecation takes place", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assemblyToDeprecate])
def cvDeprecationWriter = new ClusteredVariantDeprecationWriter(options.assemblyToDeprecate, dbEnv.mongoTemplate,
        dbEnv.submittedVariantAccessioningService,
        dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
        "EVA3196", "No SS associated with the RS's locus")
def metricCompute = dbEnv.springApplicationContext.getBean(MetricCompute.class)
[dbsnpSvoeClass, svoeClass].each{opClass ->
    new RetryableBatchingCursor<>(where("eventType").is(EventType.DEPRECATED.toString())
            .and("inactiveObjects.seq").is(options.assemblyToDeprecate)
            .and("inactiveObjects.rs").exists(true), dbEnv.mongoTemplate, opClass).each {svoes ->
        Set<ImmutablePair<String, Long>> rsHashesAndIDsToLookFor =
                svoes.collect{ EVAObjectModelUtils.toClusteredVariantEntity(
                        it.inactiveObjects[0].toSubmittedVariantEntity())}
                        .collect{new ImmutablePair<>(it.hashedMessage, it.accession)}.toSet()
        def cvesToDeprecate = [cveClass, dbsnpCveClass].collect {collectionClass ->
            dbEnv.mongoTemplate.find(query(where("asm").is(options.assemblyToDeprecate)
                    .and("_id").in(rsHashesAndIDsToLookFor.collect{it.left})
                    .and("accession").in(rsHashesAndIDsToLookFor.collect{it.right})),
                    collectionClass)
                    .findAll{rsHashesAndIDsToLookFor.contains(new ImmutablePair<>(it.hashedMessage, it.accession))}
        }.flatten()
        cvDeprecationWriter.write(cvesToDeprecate)
        metricCompute.addCount(ClusteringMetric.CLUSTERED_VARIANTS_DEPRECATED,
                cvDeprecationWriter.numDeprecatedEntities)
        metricCompute.saveMetricsCountsInDB()
}}
