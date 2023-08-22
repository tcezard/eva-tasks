package eva3363

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
import org.springframework.data.mongodb.core.query.Update
import org.springframework.transaction.PlatformTransactionManager
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
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

class DeprecatePRJEB55361IncorrectContig {

    static void deprecateSS(EVADatabaseEnvironment dbEnv, Collection<? extends SubmittedVariantEntity> svesToDeprecate) {
        if (svesToDeprecate.size() == 0) return
        def inputParameters = dbEnv.springApplicationContext.getBean(InputParameters.class)
        def svDeprecationWriter = new SubmittedVariantDeprecationWriter(inputParameters.assemblyAccession,
                dbEnv.mongoTemplate,
                dbEnv.submittedVariantAccessioningService, dbEnv.clusteredVariantAccessioningService,
                dbEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class),
                dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
                "EVA3363", "Variant remapped to a non-existent contig in the assembly")
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

def impactedAssembly = "GCA_000247795.2"
def sourceAssembly = "GCA_002263795.2"
def impactedContig = "CM008173.2"

def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class,
        ["parameters.assemblyAccession": impactedAssembly])
def impactedVariants = [sveClass, dbsnpSveClass].collect {
    new RetryableBatchingCursor(where("seq").is(impactedAssembly)
            .and("contig").is(impactedContig).and("remappedFrom").is(sourceAssembly), prodEnv.mongoTemplate, it)
}

impactedVariants.each{it.each{sves ->
    // Deprecate remapped SS with the impacted contig - see EVA-3363
    DeprecatePRJEB55361IncorrectContig.deprecateSS(prodEnv, sves)
    def impactedAccessions = sves.collect{it.accession}.toSet()
    // Remove any previously back-propagated RS to the source SS
    def backPropOps = [svoeClass, dbsnpSvoeClass].collect{collectionClass ->
        prodEnv.mongoTemplate.find(query(where("inactiveObjects.seq").is(sourceAssembly)
                .and("eventType").is(EventType.RS_BACK_PROPAGATED.toString())
                .and("accession").in(impactedAccessions)
    ), collectionClass)}.flatten()
    def ssHashesThatNeedUpdates = backPropOps.collect{it.inactiveObjects[0].hashedMessage}.toSet()
    [sveClass, dbsnpSveClass].each{prodEnv.mongoTemplate.updateMulti(collectionClass ->
            query(where("_id").in(ssHashesThatNeedUpdates)), Update.unset("backPropRS"), collectionClass)}
    ssHashesThatNeedUpdates.each{println("Removed backPropRS for SS hash $it...")}
}}
