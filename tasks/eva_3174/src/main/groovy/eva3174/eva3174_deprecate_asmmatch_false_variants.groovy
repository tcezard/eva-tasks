package eva3174

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
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script deprecates variants having asmMatch false in dbsnpSVE collection and also the corresponding RS (if any)
def cli = new CliBuilder()
cli.propertiesFile(args: 1, "Properties file to use for deprecation", required: true)
cli.assemblyToDeprecate(args: 1, "Assembly where variants with asmMatch=false should be deprecated", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

class DeprecateAsmMatchFalseSS {

    static void deprecateSS(EVADatabaseEnvironment dbEnv, Collection<? extends SubmittedVariantEntity> svesToDeprecate) {
        if (svesToDeprecate.size() == 0) return
        def inputParameters = dbEnv.springApplicationContext.getBean(InputParameters.class)
        def svDeprecationWriter = new SubmittedVariantDeprecationWriter(inputParameters.assemblyAccession,
                dbEnv.mongoTemplate,
                dbEnv.submittedVariantAccessioningService, dbEnv.clusteredVariantAccessioningService,
                dbEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class),
                dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
                "EVA3174", "Variant deprecated due to asmMatch=false")
        def dbEnvProgressListener =
                new DeprecationStepProgressListener(svDeprecationWriter, dbEnv.springApplicationContext.getBean(MetricCompute.class))

        def dbEnvJobRepository = dbEnv.springApplicationContext.getBean(JobRepository.class)
        def dbEnvTxnMgr = dbEnv.springApplicationContext.getBean(PlatformTransactionManager.class)

        def dbEnvJobBuilderFactory = new JobBuilderFactory(dbEnvJobRepository)
        def dbEnvStepBuilderFactory = new StepBuilderFactory(dbEnvJobRepository, dbEnvTxnMgr)

        def mapWtSSIterator = svesToDeprecate.iterator()
        def svDeprecateJobSteps = dbEnvStepBuilderFactory.get("stepsForSSDeprecation").chunk(new SimpleCompletionPolicy(inputParameters.chunkSize)).reader(
                new ItemStreamReader<SubmittedVariantEntity>() {
                    @Override
                    SubmittedVariantEntity read() throws Exception, UnexpectedInputException, ParseException, NonTransientResourceException {
                        return mapWtSSIterator.hasNext() ? mapWtSSIterator.next() : null
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

def dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assemblyToDeprecate])
def asmMatchFalseVariantCursor = [sveClass, dbsnpSveClass].collect{new RetryableBatchingCursor(where("seq")
        .is(options.assemblyToDeprecate).and("asmMatch").is(false), dbEnv.mongoTemplate, it)}
asmMatchFalseVariantCursor.each{it.each{sves -> DeprecateAsmMatchFalseSS.deprecateSS(dbEnv, sves)}}
