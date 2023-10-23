package eva3417


import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory
import org.springframework.batch.core.JobParameter
import org.springframework.batch.core.JobParameters
import org.springframework.batch.core.StepExecutionListener
import org.springframework.batch.core.configuration.annotation.JobBuilderFactory
import org.springframework.batch.core.configuration.annotation.StepBuilderFactory
import org.springframework.batch.core.launch.JobLauncher
import org.springframework.batch.core.launch.support.RunIdIncrementer
import org.springframework.batch.core.repository.JobRepository
import org.springframework.batch.item.*
import org.springframework.batch.repeat.policy.SimpleCompletionPolicy
import org.springframework.data.mongodb.core.BulkOperations
import org.springframework.data.mongodb.core.query.Criteria
import org.springframework.data.mongodb.core.query.Query
import org.springframework.data.mongodb.core.query.Update
import org.springframework.transaction.PlatformTransactionManager
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.batch.io.SubmittedVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.summary.SubmittedVariantSummaryFunction
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.accession.deprecate.batch.listeners.DeprecationStepProgressListener
import uk.ac.ebi.eva.accession.deprecate.parameters.InputParameters
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.metrics.metric.MetricCompute

import java.text.ParseException
import java.time.LocalDateTime
import java.util.function.Function
import java.util.stream.Collectors

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

def cli = new CliBuilder()
cli.prodPropertiesFile(args: 1, "Production properties file for connecting to prod db", required: true)
cli.hashcollisionFile(args: 2, "Path to the file to store hash collision ids", required: true)

def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

// call to remediate lowercase nucleotide
new RemediateLowerCaseNucleotide(options.prodPropertiesFile, options.hashcollisionFile).remediate()

class RemediateLowerCaseNucleotide {
    private String hashCollisionFilePath
    private EVADatabaseEnvironment dbEnv
    private Function<ISubmittedVariant, String> hashingFunction

    static def logger = LoggerFactory.getLogger(Application.class)

    RemediateLowerCaseNucleotide(prodPropertiesFile, hashcollisionFile) {
        this.dbEnv = createFromSpringContext(prodPropertiesFile, Application.class)
        this.hashCollisionFilePath = hashcollisionFile
        this.hashingFunction = new SubmittedVariantSummaryFunction().andThen(new SHA1HashingFunction())
    }

    List<SubmittedVariantEntity> getImpactedSVEList(List<SubmittedVariantEntity> sveList) {
        List<SubmittedVariantEntity> impactedSVEList = new ArrayList<>()
        for (SubmittedVariantEntity sve : sveList) {
            String ref = sve.getReferenceAllele()
            String alt = sve.getAlternateAllele()
            if (ref != ref.toUpperCase() || alt != alt.toUpperCase()) {
                impactedSVEList.add(sve)
            }
        }

        return impactedSVEList
    }


    String getSVENewHash(SubmittedVariantEntity sve) {
        SubmittedVariant temp = sve.getModel()
        temp.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
        temp.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
        return hashingFunction.apply(temp)
    }

    List<String> getExistingHashList(Collection<String> newHashList) {
        List<SubmittedVariantEntity> existingHashSVEList = dbEnv.mongoTemplate
                .find(query(where("_id").in(newHashList)), sveClass)

        return existingHashSVEList.stream()
                .map(sve -> sve.getHashedMessage())
                .collect(Collectors.toList())
    }

    void insertSVEWithNewHash(List<SubmittedVariantEntity> sveList, Map<String, String> mapOldHashNewHash) {
        List<SubmittedVariantEntity> sveWithNewHashList = new ArrayList<>()
        for (SubmittedVariantEntity sve : sveList) {
            // create a new SVE with updated ref, alt and hash
            SubmittedVariant sveModel = sve.getModel()
            sveModel.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
            sveModel.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
            SubmittedVariantEntity submittedVariantEntity = new SubmittedVariantEntity(sve.getAccession(),
                    mapOldHashNewHash.get(sve.getHashedMessage()), sveModel, sve.getVersion())

            sveWithNewHashList.add(submittedVariantEntity)
        }

        // insert updated SVEs into db
        dbEnv.mongoTemplate.insert(sveWithNewHashList, sveClass)
    }

    void deprecateSVEWithOldHash(List<SubmittedVariantEntity> sveList) {
        // create a map of sve with key being the assembly
        Map<String, List<SubmittedVariantEntity>> mapAssemblySVE = sveList.stream()
                .collect(Collectors.groupingBy(sve -> sve.getReferenceSequenceAccession()))

        // call the deprecationWriter for each assembly once
        for (Map.Entry<String, List<SubmittedVariantEntity>> entry : mapAssemblySVE.entrySet()) {
            //assembly for deprecation
            String assemblyAccession = entry.getKey()
            // SVE in assembly to deprecate
            List<SubmittedVariantEntity> svesToDeprecate = entry.getValue()

            def inputParameters = dbEnv.springApplicationContext.getBean(InputParameters.class)
            def svDeprecationWriter = new SubmittedVariantDeprecationWriter(assemblyAccession,
                    dbEnv.mongoTemplate,
                    dbEnv.submittedVariantAccessioningService, dbEnv.clusteredVariantAccessioningService,
                    dbEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class),
                    dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
                    "EVA3417", "Deprecated as nucleotide was lower case")
            def dbEnvProgressListener =
                    new DeprecationStepProgressListener(svDeprecationWriter, dbEnv.springApplicationContext.getBean(MetricCompute.class))

            def dbEnvJobRepository = dbEnv.springApplicationContext.getBean(JobRepository.class)
            def dbEnvTxnMgr = dbEnv.springApplicationContext.getBean(PlatformTransactionManager.class)

            def dbEnvJobBuilderFactory = new JobBuilderFactory(dbEnvJobRepository)
            def dbEnvStepBuilderFactory = new StepBuilderFactory(dbEnvJobRepository, dbEnvTxnMgr)

            def mapWtSSIterator = svesToDeprecate.iterator()
            def svDeprecateJobSteps = dbEnvStepBuilderFactory.get("stepsForSSDeprecation")
                    .chunk(new SimpleCompletionPolicy(inputParameters.chunkSize)).reader(
                    new ItemStreamReader<SubmittedVariantEntity>() {
                        @Override
                        SubmittedVariantEntity read() throws Exception, UnexpectedInputException, ParseException, NonTransientResourceException {
                            return mapWtSSIterator.hasNext() ? mapWtSSIterator.next() : null
                        }

                        @Override
                        void open(ExecutionContext executionContext) throws ItemStreamException {}

                        @Override
                        void update(ExecutionContext executionContext) throws ItemStreamException {}

                        @Override
                        void close() throws ItemStreamException {}
                    })
                    .writer(svDeprecationWriter)
                    .listener((StepExecutionListener) dbEnvProgressListener)
                    .build()

            def svDeprecationJob = dbEnvJobBuilderFactory.get("deprecationJob")
                    .start(svDeprecateJobSteps)
                    .incrementer(new RunIdIncrementer()).build()
            def dbEnvJobLauncher = dbEnv.springApplicationContext.getBean(JobLauncher.class)
            dbEnvJobLauncher.run(svDeprecationJob, new JobParameters(["executionDate": new JobParameter(LocalDateTime.now().toDate())]))
        }
    }

    void updateFileWithHashCollisionList(List<String> sveHashList) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(this.hashCollisionFilePath, true))) {
            for (String line : sveHashList) {
                writer.write(line);
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    void updateSVOE(List<SubmittedVariantEntity> impactedSVEList, Map<String, String> mapOldHashNewHash) {
        List<String> impactedSVEIds = impactedSVEList.stream()
                .map { sve -> sve.getHashedMessage() }
                .collect(Collectors.toList())
        Criteria filterCriteria = new Criteria("inactiveObjects.hashedMessage").in(impactedSVEIds)
                .and("eventType").is("UPDATED")
        new RetryableBatchingCursor<SubmittedVariantOperationEntity>(filterCriteria, dbEnv.mongoTemplate, svoeClass).each {
            sveOperationBatchList ->
                {
                    def svoeBulkUpdates = dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, svoeClass)

                    for (SubmittedVariantOperationEntity svoe : sveOperationBatchList) {
                        // create update statement
                        SubmittedVariantInactiveEntity inactiveSVEDocument = svoe.getInactiveObjects()[0]
                        String oldHash = inactiveSVEDocument.getHashedMessage()
                        String newRef = inactiveSVEDocument.getReferenceAllele().toUpperCase()
                        String newAlt = inactiveSVEDocument.getAlternateAllele().toUpperCase()
                        String newHash = mapOldHashNewHash.get(oldHash)

                        Query findQuery = query(where("inactiveObjects.hashedMessage").is(oldHash))
                        Update updateQuery = new Update()
                                .set("inactiveObjects.\$.ref", newRef)
                                .set("inactiveObjects.\$.alt", newAlt)
                                .set("inactiveObjects.\$.hashedMessage", newHash)

                        svoeBulkUpdates.updateOne(findQuery, updateQuery)
                    }

                    //execute all bulk updates
                    svoeBulkUpdates.execute()
                }
        }

    }

    void remediate() {
        // read all the documents from SubmittedVariantEntity
        new RetryableBatchingCursor<SubmittedVariantEntity>(new Criteria(), dbEnv.mongoTemplate, sveClass).each {
            sveBatchList ->
                {
                    // find which sve from the batch are impacted by lower case nucleotide issue
                    List<SubmittedVariantEntity> impactedSVEList = getImpactedSVEList(sveBatchList)

                    // if there are any sve impacted
                    if (!impactedSVEList.isEmpty()) {
                        // make a map of old hash ids of sve and new hash
                        Map<String, String> mapOldHashNewHash = impactedSVEList.stream()
                                .collect(Collectors.toMap(sve -> sve.getHashedMessage(), sve -> getSVENewHash(sve)))

                        // find which of the new hashes already exist in db (will cause hash collision)
                        List<String> newHashExistList = getExistingHashList(mapOldHashNewHash.values())

                        // partition the impacted sve list into 2 parts  (hash collision and no hash collision)
                        Map<Boolean, List<SubmittedVariantEntity>> svePartitionMap = impactedSVEList.stream()
                                .collect(Collectors.partitioningBy(sve -> newHashExistList.contains(mapOldHashNewHash.get(sve.getHashedMessage()))))
                        List<SubmittedVariantEntity> noHashCollisionSVEList = svePartitionMap.get(Boolean.FALSE)

                        // sve with no hash collision - insert sve with new hash and deprecate sve with old hash
                        // can't update in place as _id is immutable
                        logger.info("Impacted sve List (No Hash Collision): " + noHashCollisionSVEList)
                        insertSVEWithNewHash(noHashCollisionSVEList, mapOldHashNewHash)
                        deprecateSVEWithOldHash(noHashCollisionSVEList)

                        List<SubmittedVariantEntity> hashCollisionSVEList = svePartitionMap.get(Boolean.TRUE)
                        if (!hashCollisionSVEList.isEmpty()) {
                            logger.info("Impacted sve List (Hash Collision): " + hashCollisionSVEList)
                            List<String> collisionList = hashCollisionSVEList.stream()
                                    .map(sve -> sve.getHashedMessage())
                                    .collect(Collectors.toList())

                            updateFileWithHashCollisionList(collisionList)
                            deprecateSVEWithOldHash(hashCollisionSVEList)
                        }

                        //update SVE Operations
                        updateSVOE(impactedSVEList, mapOldHashNewHash)

                    }
                }
        }
    }
}

