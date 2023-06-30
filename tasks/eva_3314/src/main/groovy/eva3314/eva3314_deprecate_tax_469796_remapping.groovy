package eva3309

import groovy.cli.picocli.CliBuilder
import org.springframework.batch.core.JobParameter
import org.springframework.batch.core.JobParameters
import org.springframework.batch.core.StepExecutionListener
import org.springframework.batch.core.configuration.annotation.JobBuilderFactory
import org.springframework.batch.core.configuration.annotation.StepBuilderFactory
import org.springframework.batch.core.launch.JobLauncher
import org.springframework.batch.core.launch.support.RunIdIncrementer
import org.springframework.batch.core.repository.JobRepository
import org.springframework.batch.item.ExecutionContext
import org.springframework.batch.item.ItemStreamException
import org.springframework.batch.item.ItemStreamReader
import org.springframework.batch.item.NonTransientResourceException
import org.springframework.batch.item.UnexpectedInputException
import org.springframework.batch.repeat.policy.SimpleCompletionPolicy
import org.springframework.transaction.PlatformTransactionManager

import uk.ac.ebi.eva.accession.core.batch.io.SubmittedVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.accession.deprecate.batch.listeners.DeprecationStepProgressListener
import uk.ac.ebi.eva.accession.deprecate.parameters.InputParameters
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.metrics.metric.MetricCompute

import java.text.ParseException
import java.time.LocalDateTime

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

class Deprecate469796TaxRemappedSS {

    static void deprecateSS(EVADatabaseEnvironment dbEnv, Collection<? extends SubmittedVariantEntity> svesToDeprecate) {
        if (svesToDeprecate.size() == 0) return
        def inputParameters = dbEnv.springApplicationContext.getBean(InputParameters.class)
        def svDeprecationWriter = new SubmittedVariantDeprecationWriter(inputParameters.assemblyAccession,
                dbEnv.mongoTemplate,
                dbEnv.submittedVariantAccessioningService, dbEnv.clusteredVariantAccessioningService,
                dbEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class),
                dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
                "EVA3314", "Variant incorrectly remapped to the wrong assembly for tax 469796")
        def dbEnvProgressListener =
                new DeprecationStepProgressListener(svDeprecationWriter, dbEnv.springApplicationContext.getBean(MetricCompute.class))

        def dbEnvJobRepository = dbEnv.springApplicationContext.getBean(JobRepository.class)
        def dbEnvTxnMgr = dbEnv.springApplicationContext.getBean(PlatformTransactionManager.class)

        def dbEnvJobBuilderFactory = new JobBuilderFactory(dbEnvJobRepository)
        def dbEnvStepBuilderFactory = new StepBuilderFactory(dbEnvJobRepository, dbEnvTxnMgr)

        def ssIterator = svesToDeprecate.iterator()
        def svDeprecateJobSteps = dbEnvStepBuilderFactory.get("stepsForSSDeprecation").chunk(new SimpleCompletionPolicy(inputParameters.chunkSize)).reader(
                new ItemStreamReader<SubmittedVariantEntity>() {
                    @Override
                    SubmittedVariantEntity read() throws Exception, UnexpectedInputException, ParseException, NonTransientResourceException {
                        return ssIterator.hasNext() ? ssIterator.next() : null
                    }

                    @Override
                    void open(ExecutionContext executionContext) throws ItemStreamException {

                    }

                    @Override
                    void update(ExecutionContext executionContext) throws ItemStreamException {

                    }

                    @Override
                    void close() throws ItemStreamException {

                    }
                }).writer(svDeprecationWriter).listener((StepExecutionListener) dbEnvProgressListener).build()

        def svDeprecationJob = dbEnvJobBuilderFactory.get("deprecationJob").start(svDeprecateJobSteps).incrementer(
                new RunIdIncrementer()).build()
        def dbEnvJobLauncher = dbEnv.springApplicationContext.getBean(JobLauncher.class)
        dbEnvJobLauncher.run(svDeprecationJob, new JobParameters(["executionDate": new JobParameter(
                LocalDateTime.now().toDate())]))
    }
}

// This script deprecates variants incorrectly remapped in the Ovis orientalis species (tax 469796)
def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file to use for deprecation", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def impactedAssembly = "GCA_016772045.1"
def impactedTax = 469796

def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class,
        ["parameters.assemblyAccession": impactedAssembly])
def tax469796VariantCursor = [sveClass, dbsnpSveClass].collect {
    new RetryableBatchingCursor(where("seq").is(impactedAssembly)
            .and("tax").is(impactedTax).and("remappedFrom").exists(true), prodEnv.mongoTemplate, it)
}
tax469796VariantCursor.each{it.each{sves ->
    // Deprecate matching SS with no duplicates in the remapped assembly
    Deprecate469796TaxRemappedSS.deprecateSS(prodEnv, sves)
}}
// Delete all operations involving the taxonomy in GCA_016772045.1
def tax469796SubmittedOpCursor = [svoeClass, dbsnpSvoeClass].collect {
    new RetryableBatchingCursor<>(where("inactiveObjects.seq").is(impactedAssembly)
            .and("inactiveObjects.tax").is(impactedTax), prodEnv.mongoTemplate, it)
}
int numSubmittedOpsRemoved = 0
tax469796SubmittedOpCursor.each {currCursor -> currCursor.each{svoes ->
    def impactedOpIds = svoes.collect{it.id}.toSet()
    def deletedRecords = prodEnv.mongoTemplate.findAllAndRemove(query(where("_id").in(impactedOpIds)),
            currCursor.collectionClass)
    numSubmittedOpsRemoved += deletedRecords.size()
    println("${numSubmittedOpsRemoved} operations removed so far " +
            "in collection ${currCursor.collectionClass.simpleName}...")
}}
