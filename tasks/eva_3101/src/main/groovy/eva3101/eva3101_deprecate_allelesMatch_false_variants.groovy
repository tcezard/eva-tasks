package eva3101

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
import org.springframework.batch.item.ParseException
import org.springframework.batch.item.UnexpectedInputException
import org.springframework.batch.repeat.policy.SimpleCompletionPolicy
import org.springframework.transaction.PlatformTransactionManager
import uk.ac.ebi.eva.accession.core.batch.io.SubmittedVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.accession.deprecate.batch.listeners.DeprecationStepProgressListener
import uk.ac.ebi.eva.accession.deprecate.parameters.InputParameters
import uk.ac.ebi.eva.groovy.commons.RetryableCursor
import uk.ac.ebi.eva.metrics.metric.MetricCompute

import java.time.LocalDateTime

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

// This script populates the development environment with variants having allelesMatch false in dbsnpSVE collection and also the corresponding RS (if any)
def cli = new CliBuilder()
cli.propertiesFile(args: 1, "Properties file for accessioning", required: true)
cli.assemblyToDeprecate(args: 1, "Assembly where variants with allelesMatch=false should be deprecated", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["parameters.assemblyAccession": options.assemblyToDeprecate])

def dbEnvDeprecationWriter = new SubmittedVariantDeprecationWriter(options.assemblyToDeprecate, dbEnv.mongoTemplate,
        dbEnv.submittedVariantAccessioningService, dbEnv.clusteredVariantAccessioningService,
        dbEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class),
        dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
        "EVA3101", "Variant deprecated due to alleles mismatch with the reference assembly.")
def dbEnvProgressListener =
        new DeprecationStepProgressListener(dbEnvDeprecationWriter, dbEnv.springApplicationContext.getBean(MetricCompute.class))

def inputParameters = dbEnv.springApplicationContext.getBean(InputParameters.class)
def dbEnvJobRepository = dbEnv.springApplicationContext.getBean(JobRepository.class)
def dbEnvTxnMgr = dbEnv.springApplicationContext.getBean(PlatformTransactionManager.class)

def dbEnvJobBuilderFactory = new JobBuilderFactory(dbEnvJobRepository)
def dbEnvStepBuilderFactory = new StepBuilderFactory(dbEnvJobRepository, dbEnvTxnMgr)

def allelesMatchFalseVariantIterator = new RetryableCursor(where("seq")
        .is(options.assemblyToDeprecate).and("allelesMatch").exists(true),
        dbEnv.mongoTemplate, dbsnpSveClass).iterator()
def dbEnvJobSteps = dbEnvStepBuilderFactory.get("stepsForEVA3101Deprecation").chunk(
        new SimpleCompletionPolicy(inputParameters.chunkSize)).reader(
        new ItemStreamReader<SubmittedVariantEntity>() {
            @Override
            SubmittedVariantEntity read() throws Exception, UnexpectedInputException, ParseException, NonTransientResourceException {
                return allelesMatchFalseVariantIterator.hasNext()? allelesMatchFalseVariantIterator.next(): null
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
        }).writer(dbEnvDeprecationWriter).listener((StepExecutionListener)dbEnvProgressListener).build()


def dbEnvDeprecationJob = dbEnvJobBuilderFactory.get("deprecationJobSteps").start(dbEnvJobSteps).incrementer(
        new RunIdIncrementer()).build()
def dbEnvJobLauncher = dbEnv.springApplicationContext.getBean(JobLauncher.class)
dbEnvJobLauncher.run(dbEnvDeprecationJob, new JobParameters(["executionDate": new JobParameter(
        LocalDateTime.now().toDate())]))
