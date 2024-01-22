package eva3417

import com.google.gson.Gson
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
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.batch.io.SubmittedVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantEntity
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

import java.time.LocalDateTime
import java.util.function.Function
import java.util.stream.Collectors

import static org.junit.Assert.assertEquals
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

def cli = new CliBuilder()
cli.envPropertiesFile(args: 1, "properties file for connecting to db", required: true)
cli.summaryFilePath(args: 2, "Path to the summary file containing the involved sve details", required: true)

def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

new BothSVEHasValidRSRemediation(options.envPropertiesFile, options.summaryFilePath).remediate()

class BothSVEHasValidRSRemediation {
    private EVADatabaseEnvironment dbEnv
    private String summaryFilePath
    private Function<ISubmittedVariant, String> hashingFunction
    static def logger = LoggerFactory.getLogger(BothSVEHasValidRSRemediation.class)
    private Gson gson

    BothSVEHasValidRSRemediation(envPropertiesFile, summaryFilePath) {
        this.dbEnv = createFromSpringContext(envPropertiesFile, Application.class)
        this.summaryFilePath = summaryFilePath
        this.hashingFunction = new SubmittedVariantSummaryFunction().andThen(new SHA1HashingFunction())
        this.gson = new Gson()
    }

    void remediate() {
        // read collision details from summary file
        List<BothRSValidCollision> bothRSValidCollisionList
        try (BufferedReader br = new BufferedReader(new FileReader(this.summaryFilePath))) {
            List<String> lines = br.readLines()
            // read only those collisions where both sve has valid rs
            bothRSValidCollisionList = lines.stream()
                    .skip(1)
                    .map(line -> new BothRSValidCollision(line))
                    .filter(cd -> cd.getCategory() == "BOTH_SVE_HAS_RS_BOTH_VALID")
        }

        // get involved SVEs and CVEs
        Map<String, SubmittedVariantEntity> mapOfHashAndSVE = getMapOfHashAndSVE(bothRSValidCollisionList)
        Map<Long, ClusteredVariantEntity> mapOfAccessionAndCVE = getMapOfAccessionAndCVE(bothRSValidCollisionList)

        // Remediate the case where both SVE has valid RS.
        // Even though both the RS are valid/exist, they have different start and only one of start matches with start in SVE
        // Keep the SVE with RS whose start matches with start of the SVE (remediate if necessary)
        // deprecate the SVE with RS whose start does not match with start of the SVE
        List<SubmittedVariantEntity> sveInsertList = new ArrayList<>()
        List<SubmittedVariantEntity> sveDeprecateList = new ArrayList<>()
        for (BothRSValidCollision collision : bothRSValidCollisionList) {
            SubmittedVariantEntity fileSVE = mapOfHashAndSVE.get(collision.getSveInFile())
            SubmittedVariantEntity dbSVE = mapOfHashAndSVE.get(collision.getSveInDB())
            ClusteredVariantEntity fileCVE = mapOfAccessionAndCVE.get(fileSVE.getClusteredVariantAccession())
            ClusteredVariantEntity dbCVE = mapOfAccessionAndCVE.get(dbSVE.getClusteredVariantAccession())

            // both SVE has same start, we can take it from any
            Long sveStart = fileSVE.getStart()
            if (sveStart == fileCVE.getStart() && sveStart == dbCVE.getStart()) {
                // Could not find this case in investigation, but logging this check just in case
                logger.info("SVE_START_MATCHES_WITH_BOTH_RS")
                continue
            } else if (sveStart == fileCVE.getStart() || sveStart == dbCVE.getStart()) {
                logger.info(sveStart == fileCVE.getStart() ? "SVE_START_MATCHES_WITH_FILE_RS" : "SVE_START_MATCHES_WITH_DB_RS")
                // assume fileSVE has the RS with correct start and needs to be kept
                SubmittedVariantEntity sveToKeep = fileSVE
                SubmittedVariantEntity sveToDeprecate = dbSVE
                if (sveStart == dbCVE.getStart()) {
                    sveToKeep = dbSVE
                    sveToDeprecate = fileSVE
                }
                sveInsertList.add(sveToKeep)
                sveDeprecateList.add(sveToDeprecate)
            } else {
                // Could not find this case in investigation, but logging this check just in case
                logger.info("SVE_START_MATCHES_WITH_NONE_RS")
                continue
            }
        }

        // deprecate SVE
        logger.info("List of SVE to deprecate: " + gson.toJson(sveDeprecateList))
        deprecateSVE(sveDeprecateList)

        // insert SVE with new hash (after remediation)
        insertSVEWithNewHash(sveInsertList, dbsnpSveClass)
        // insert update operatons (for remediation done)
        insertSVOEUpdateOp(sveInsertList, dbsnpSvoeClass)
        // remove existing sve (with lowercase nucleotide)
        logger.info("List of SVE to remove: " + gson.toJson(sveInsertList))
        removeImpactedSVE(sveInsertList, dbsnpSveClass)

        // update existing SVOE where a lowercase nucleotide was involved
        Map<String, String> mapOldHashNewHash = sveInsertList.stream()
                .collect(Collectors.toMap(sve -> sve.getHashedMessage(), sve -> getSVENewHash(sve)))
        updateExistingSVOE(sveInsertList, mapOldHashNewHash, dbsnpSvoeClass)
    }

    void deprecateSVE(List<SubmittedVariantEntity> sveList) {
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
                    "EVA3472", "Deprecated as collision with lowercase nucleotide - RS associated has wrong start")
            def dbEnvProgressListener =
                    new DeprecationStepProgressListener(svDeprecationWriter, dbEnv.springApplicationContext.getBean(MetricCompute.class))

            def dbEnvJobRepository = dbEnv.springApplicationContext.getBean(JobRepository.class)
            def dbEnvTxnMgr = dbEnv.springApplicationContext.getBean(PlatformTransactionManager.class)

            def dbEnvJobBuilderFactory = new JobBuilderFactory(dbEnvJobRepository)
            def dbEnvStepBuilderFactory = new StepBuilderFactory(dbEnvJobRepository, dbEnvTxnMgr)

            def mapWtSSIterator = svesToDeprecate.iterator()
            def svDeprecateJobSteps = dbEnvStepBuilderFactory.get("stepsForSSDeprecation")
                    .chunk(new SimpleCompletionPolicy(inputParameters.chunkSize)).reader(new ItemStreamReader<SubmittedVariantEntity>() {
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

    void insertSVEWithNewHash(List<SubmittedVariantEntity> sveList, Class ssClass) {
        List<SubmittedVariantEntity> sveWithNewHashList = new ArrayList<>()
        for (SubmittedVariantEntity sve : sveList) {
            // create a new SVE with updated ref, alt and hash
            SubmittedVariant sveModel = sve.getModel()
            sveModel.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
            sveModel.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
            SubmittedVariantEntity submittedVariantEntity = new SubmittedVariantEntity(sve.getAccession(),
                    getSVENewHash(sve), sveModel, sve.getVersion())

            sveWithNewHashList.add(submittedVariantEntity)
        }

        // insert updated SVEs into db
        logger.info("List of SVE to insert: " + gson.toJson(sveWithNewHashList))
        dbEnv.mongoTemplate.insert(sveWithNewHashList, ssClass)
    }

    String getSVENewHash(SubmittedVariantEntity sve) {
        SubmittedVariant temp = sve.getModel()
        temp.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
        temp.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
        return hashingFunction.apply(temp)
    }

    void insertSVOEUpdateOp(List<SubmittedVariantEntity> sveList, Class ssClass) {
        def sveUpdateOpList = new ArrayList<>()
        for (SubmittedVariantEntity sve : sveList) {
            def updateOpToWrite = ssClass.equals(sveClass) ? new SubmittedVariantOperationEntity() : new DbsnpSubmittedVariantOperationEntity()
            updateOpToWrite.fill(EventType.UPDATED, sve.getAccession(), null,
                    "EVA3472 - SS updated with upper case ref, alt and a new hash",
                    Arrays<SubmittedVariantInactiveEntity>.asList(new SubmittedVariantInactiveEntity(sve)).asList())

            updateOpToWrite.setId("EVA3472_UPDATED_" + sve.getReferenceSequenceAccession() + "_" + sve.getAccession())

            sveUpdateOpList.add(updateOpToWrite)
        }

        logger.info("List of SVE operations to insert: " + gson.toJson(sveUpdateOpList))
        dbEnv.mongoTemplate.insert(sveUpdateOpList, ssClass.equals(sveClass) ? svoeClass : dbsnpSvoeClass)
    }

    void removeImpactedSVE(List<SubmittedVariantEntity> impactedSVEList, Class ssClass) {
        List<String> impactedSVEIds = impactedSVEList.stream()
                .map(sve -> sve.getHashedMessage())
                .collect(Collectors.toList())
        Query findQuery = query(where("_id").in(impactedSVEIds))

        dbEnv.mongoTemplate.remove(findQuery, ssClass)
    }

    void updateExistingSVOE(List<SubmittedVariantEntity> impactedSVEList, Map<String, String> mapOldHashNewHash,
                            Class ssClass) {
        List<String> impactedSVEIds = impactedSVEList.stream()
                .map(sve -> sve.getHashedMessage())
                .collect(Collectors.toList())
        List<String> impactedAssemblies = impactedSVEList.stream()
                .map(sve -> sve.getReferenceSequenceAccession())
                .collect(Collectors.toSet()).stream()
                .collect(Collectors.toList())
        List<Long> impactedAccessions = impactedSVEList.stream()
                .map(sve -> sve.getAccession())
                .collect(Collectors.toList())

        Criteria filterCriteria = new Criteria("inactiveObjects.seq").in(impactedAssemblies)
                .and("accession").in(impactedAccessions)
                .and("eventType").is("UPDATED")
                .and("inactiveObjects.hashedMessage").in(impactedSVEIds)
                .and("_id").not().regex("^EVA3472_UPDATED_")
        new RetryableBatchingCursor<SubmittedVariantOperationEntity>(filterCriteria, dbEnv.mongoTemplate, ssClass).each { sveOperationBatchList ->
            {
                logger.info("Existing SVOE operations to update: " + gson.toJson(sveOperationBatchList))

                def svoeBulkUpdates = dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, ssClass)

                for (SubmittedVariantOperationEntity svoe : sveOperationBatchList) {
                    // create update statement
                    SubmittedVariantInactiveEntity inactiveSVEDocument = svoe.getInactiveObjects()[0]
                    String oldHash = inactiveSVEDocument.getHashedMessage()
                    String newRef = inactiveSVEDocument.getReferenceAllele().toUpperCase()
                    String newAlt = inactiveSVEDocument.getAlternateAllele().toUpperCase()
                    String newHash = mapOldHashNewHash.get(oldHash)

                    Query findQuery = query(where("_id").is(svoe.getId()).and("inactiveObjects.hashedMessage").is(oldHash))
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


    Map<String, SubmittedVariantEntity> getMapOfHashAndSVE(List<BothRSValidCollision> bothRSValidCollisionsList) {
        Set<String> allSVEIds = new HashSet<>()
        Set<String> fileIds = bothRSValidCollisionsList.stream()
                .map(cd -> cd.getSveInFile()).collect(Collectors.toSet())
        Set<String> dbIds = bothRSValidCollisionsList.stream()
                .map(cd -> cd.getSveInDB()).collect(Collectors.toSet())
        allSVEIds.addAll(fileIds)
        allSVEIds.addAll(dbIds)

        List<SubmittedVariantEntity> allSVEList = getAllSVEByHash(allSVEIds)

        Map<String, SubmittedVariantEntity> mapOfHashAndSVE = allSVEList.stream()
                .collect(Collectors.toMap(sve -> sve.getHashedMessage(), sve -> sve))

        return mapOfHashAndSVE
    }

    Map<Long, ClusteredVariantEntity> getMapOfAccessionAndCVE(List<BothRSValidCollision> bothRSValidCollisionsList) {
        Set<String> assemblySet = bothRSValidCollisionsList.stream()
                .map(cd -> cd.assembly).collect(Collectors.toSet())
        Set<Long> allCVEAcc = new HashSet<>()
        Set<Long> fileAcc = bothRSValidCollisionsList.stream().map(cd -> cd.getFileRS())
                .collect(Collectors.toSet())
        Set<Long> dbAcc = bothRSValidCollisionsList.stream().map(cd -> cd.getDbRS())
                .collect(Collectors.toSet())
        allCVEAcc.addAll(fileAcc)
        allCVEAcc.addAll(dbAcc)

        List<ClusteredVariantEntity> allCVEList = getAllCVEByAccession(assemblySet, allCVEAcc)

        Map<Long, ClusteredVariantEntity> mapOfAccAndCVE = allCVEList.stream()
                .collect(Collectors.toMap(cve -> cve.getAccession(), cve -> cve))

        return mapOfAccAndCVE
    }

    List<SubmittedVariantEntity> getAllSVEByHash(Set<String> sveIdList) {
        List<SubmittedVariantEntity> sveList = [sveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(sveIdList)), ssClass))
                .flatten()
        // In investigation we found that all the SVE involved are in dbsnp, checking it here for confirmation
        assertEquals(0, sveList.size())

        List<SubmittedVariantEntity> allSVEInvolved = [dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(sveIdList)), ssClass))
                .flatten()

        return allSVEInvolved
    }

    List<ClusteredVariantEntity> getAllCVEByAccession(Set<String> assemblySet, Set<Long> cveAccessionList) {
        List<ClusteredVariantEntity> allCVEInvolved = [cveClass, dbsnpCveClass]
                .collectMany(cveClass -> dbEnv.mongoTemplate.find(query(where("asm").in(assemblySet)
                        .and("accession").in(cveAccessionList)), cveClass))
                .flatten()

        return allCVEInvolved
    }

}

class BothRSValidCollision {
    private String assembly

    private String sveInFile
    private Long fileAccession
    private Long fileRS
    private Boolean fileRSValid
    private String fileRemappedFrom

    private String sveInDB
    private Long dbAccession
    private Long dbRS
    private Boolean dbRSValid
    private String dbRemappedFrom

    private String fileRef
    private String fileAlt
    private String dbRef
    private String dbAlt

    private String category

    BothRSValidCollision(String line) {
        String[] values = line.split(",")
        if (values.length != 16) {
            throw new RuntimeException("Values are missing. Values Length: " + values.length + " Line: " + line)
        }

        this.assembly = getString(values[0])
        this.sveInFile = getString(values[1])
        this.fileAccession = getLong(values[2])
        this.fileRS = getLong(values[3])
        this.fileRSValid = Boolean.parseBoolean(values[4])
        this.fileRemappedFrom = getString(values[5])
        this.sveInDB = getString(values[6])
        this.dbAccession = getLong(values[7])
        this.dbRS = getLong(values[8])
        this.dbRSValid = Boolean.parseBoolean(values[9])
        this.dbRemappedFrom = getString(values[10])
        this.fileRef = getString(values[11])
        this.fileAlt = getString(values[12])
        this.dbRef = getString(values[13])
        this.dbAlt = getString(values[14])
        this.category = getString(values[15])
    }

    private String getString(String value) {
        if (value == null || value.isEmpty() || value == "null") {
            return ""
        } else {
            return value
        }
    }

    private Long getLong(String value) {
        if (value == null || value.isEmpty() || value == "null") {
            return null
        } else {
            return Long.parseLong(value)
        }
    }

    String getAssembly() {
        return assembly
    }

    String getSveInFile() {
        return sveInFile
    }

    Long getFileAccession() {
        return fileAccession
    }

    Long getFileRS() {
        return fileRS
    }

    Boolean getFileRSValid() {
        return fileRSValid
    }

    String getFileRemappedFrom() {
        return fileRemappedFrom
    }

    String getSveInDB() {
        return sveInDB
    }

    Long getDbAccession() {
        return dbAccession
    }

    Long getDbRS() {
        return dbRS
    }

    Boolean getDbRSValid() {
        return dbRSValid
    }

    String getDbRemappedFrom() {
        return dbRemappedFrom
    }

    String getFileRef() {
        return fileRef
    }

    String getFileAlt() {
        return fileAlt
    }

    String getDbRef() {
        return dbRef
    }

    String getDbAlt() {
        return dbAlt
    }

    String getCategory() {
        return category
    }
}