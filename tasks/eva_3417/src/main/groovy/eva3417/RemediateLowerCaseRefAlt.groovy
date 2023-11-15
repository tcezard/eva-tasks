package eva3417

import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.BulkOperations
import org.springframework.data.mongodb.core.query.Criteria
import org.springframework.data.mongodb.core.query.Query
import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.summary.SubmittedVariantSummaryFunction
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import java.util.function.Function
import java.util.stream.Collectors

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

def cli = new CliBuilder()
cli.envPropertiesFile(args: 1, "Production properties file for connecting to prod db", required: true)
cli.hashcollisionFile(args: 2, "Path to the file to store hash collision ids", required: true)

def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

// call to remediate lowercase nucleotide
new RemediateLowerCaseNucleotide(options.envPropertiesFile, options.hashcollisionFile).remediate()

class RemediateLowerCaseNucleotide {
    private String hashCollisionFilePath
    private EVADatabaseEnvironment dbEnv
    private Function<ISubmittedVariant, String> hashingFunction
    static def logger = LoggerFactory.getLogger(Application.class)

    RemediateLowerCaseNucleotide(envPropertiesFile, hashcollisionFile) {
        this.dbEnv = createFromSpringContext(envPropertiesFile, Application.class)
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

    Set<SubmittedVariantEntity> getExistingHashSVEList(Collection<String> newHashList) {
        List<SubmittedVariantEntity> existingHashSVEList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(newHashList)), ssClass))
                .flatten()

        return existingHashSVEList
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

    void updateFileWithHashCollisionList(List<SubmittedVariantEntity> sveList, Map<String, String> mapOldHashNewHash,
                                         Set<SubmittedVariantEntity> alreadyExistingHashSVEList) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(this.hashCollisionFilePath, true))) {
            for (SubmittedVariantEntity sve : sveList) {
                String sveNewHash = mapOldHashNewHash.get(sve.getHashedMessage())
                SubmittedVariantEntity existingSVE = alreadyExistingHashSVEList.stream()
                        .filter(s -> s.getHashedMessage().equals(sveNewHash))
                        .findFirst().get()
                //capture assembly, hash and accession of both the sve that has a collision
                String line = sve.getReferenceSequenceAccession() + "," + sve.getHashedMessage() + "," + sve.getAccession() + "," + existingSVE.getAccession()
                writer.write(line)
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    void updateExistingSVOE(List<SubmittedVariantEntity> impactedSVEList, Map<String, String> mapOldHashNewHash) {
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


    void insertSVOEUpdateOp(List<SubmittedVariantEntity> sveList, Class ssClass) {
        def sveUpdateOpList = new ArrayList<>()
        for (SubmittedVariantEntity sve : sveList) {
            def updateOpToWrite = ssClass.equals(sveClass)
                    ? new SubmittedVariantOperationEntity() : new DbsnpSubmittedVariantOperationEntity()
            updateOpToWrite.fill(EventType.UPDATED, sve.getAccession(), null,
                    "EVA3417 - SS updated with upper case ref, alt and a new hash",
                    Arrays<SubmittedVariantInactiveEntity>.asList(new SubmittedVariantInactiveEntity(sve)).asList())

            updateOpToWrite.setId("EVA3417_UPDATED_${sve.getReferenceSequenceAccession()}_${sve.getHashedMessage()}")

            sveUpdateOpList.add(updateOpToWrite)
        }

        dbEnv.mongoTemplate.insert(sveUpdateOpList, ssClass.equals(sveClass) ? svoeClass : dbsnpSvoeClass)
    }


    void removeImpactedSVE(List<SubmittedVariantEntity> impactedSVEList, Class ssClass) {
        List<String> impactedSVEIds = impactedSVEList.stream()
                .map(sve -> sve.getHashedMessage())
                .collect(Collectors.toList())
        Query findQuery = query(where("_id").in(impactedSVEIds))

        dbEnv.mongoTemplate.findAllAndRemove(findQuery, ssClass)
    }


    void remediate() {
        [sveClass, dbsnpSveClass].each {
            ssClass ->
                {
                    logger.info("Started processing for: " + ssClass)
                    long totalProcessed = 0
                    long totalImpacted = 0
                    // read all the documents from SubmittedVariantEntity
                    new RetryableBatchingCursor<SubmittedVariantEntity>(new Criteria(), dbEnv.mongoTemplate, ssClass).each {
                        sveBatchList ->
                            {
                                long currentBatchSize = sveBatchList.size()
                                logger.info("processing in this batch " + currentBatchSize)

                                // find which sve from the batch are impacted by lower case nucleotide issue
                                List<SubmittedVariantEntity> impactedSVEList = getImpactedSVEList(sveBatchList)

                                // if there are any sve impacted
                                if (!impactedSVEList.isEmpty()) {
                                    logger.info("impacted in this batch = " + impactedSVEList.size())
                                    totalImpacted = totalImpacted + impactedSVEList.size()

                                    // make a map of old hash ids of sve and new hash
                                    Map<String, String> mapOldHashNewHash = impactedSVEList.stream()
                                            .collect(Collectors.toMap(sve -> sve.getHashedMessage(), sve -> getSVENewHash(sve)))

                                    // find which of the new hashes already exist in db (will cause hash collision)
                                    Set<SubmittedVariantEntity> alreadyExistingHashSVEList = getExistingHashSVEList(mapOldHashNewHash.values())
                                    Set<String> newHashExistSet = alreadyExistingHashSVEList.stream()
                                            .map(sve -> sve.getHashedMessage())
                                            .collect(Collectors.toSet())

                                    // partition the impacted sve list into 2 parts  (hash collision and no hash collision)
                                    Map<Boolean, List<SubmittedVariantEntity>> svePartitionMap = impactedSVEList.stream()
                                            .collect(Collectors.groupingBy { sve -> newHashExistSet.contains(mapOldHashNewHash.get(sve.getHashedMessage())) })

                                    // sve with no hash collision
                                    List<SubmittedVariantEntity> noHashCollisionSVEList = svePartitionMap.get(Boolean.FALSE)
                                    logger.info("Impacted sve List (No Hash Collision): " + noHashCollisionSVEList)
                                    // insert sve with new hash (can't update in place as _id is immutable)
                                    insertSVEWithNewHash(noHashCollisionSVEList, mapOldHashNewHash)
                                    // insert sve update operation
                                    insertSVOEUpdateOp(noHashCollisionSVEList, ssClass)
                                    // remove all impacted sve from the table
                                    removeImpactedSVE(noHashCollisionSVEList, ssClass)


                                    // sve with hash collision
                                    List<SubmittedVariantEntity> hashCollisionSVEList = svePartitionMap.get(Boolean.TRUE)
                                    if (hashCollisionSVEList != null && !hashCollisionSVEList.isEmpty()) {
                                        logger.info("Impacted sve List (Hash Collision): " + hashCollisionSVEList)
                                        // capture in file - sve hash collision details
                                        updateFileWithHashCollisionList(hashCollisionSVEList, mapOldHashNewHash, alreadyExistingHashSVEList)
                                    }

                                    //update existing SVE Operations
                                    updateExistingSVOE(impactedSVEList, mapOldHashNewHash)
                                }

                                totalProcessed = totalProcessed + currentBatchSize
                                logger.info("total Processed till now = " + totalProcessed)
                                logger.info("total impacted till now = " + totalImpacted)
                            }
                    }
                }
        }

    }
}

