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

import java.time.LocalDateTime
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
        Map<String, SubmittedVariantEntity> mapOfHashAndExistingSVE = alreadyExistingHashSVEList.stream()
                                    .collect(Collectors.toMap(sve->sve.getHashedMessage(), sve->sve))

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(this.hashCollisionFilePath, true))) {
            for (SubmittedVariantEntity sve : sveList) {
                String sveNewHash = mapOldHashNewHash.get(sve.getHashedMessage())
                SubmittedVariantEntity existingSVE = mapOfHashAndExistingSVE.get(sveNewHash)
                //capture assembly, hash and accession of both the sve that has a collision
                String line = sve.getReferenceSequenceAccession() + "," + sve.getHashedMessage() + "," + existingSVE.getHashedMessage() + "," + sve.getAccession() + "," + existingSVE.getAccession()
                writer.write(line)
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    void updateExistingSVOE(List<SubmittedVariantEntity> impactedSVEList, Map<String, String> mapOldHashNewHash) {
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
                .and("_id").not().regex("^EVA3417_UPDATED_");
        new RetryableBatchingCursor<SubmittedVariantOperationEntity>(filterCriteria, dbEnv.mongoTemplate, svoeClass).each {
            sveOperationBatchList ->
                {
                    logger.info("Existing SVOE operations to update: " + sveOperationBatchList)

                    def svoeBulkUpdates = dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, svoeClass)

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

        dbEnv.mongoTemplate.remove(findQuery, ssClass)
    }


    List<SubmittedVariantEntity> findHashCollisionInBatch(List<SubmittedVariantEntity> sveList,
                                                          Map<String,String> mapOldHashNewHash){
        List<SubmittedVariantEntity> hashCollisionSVEInBatch = new ArrayList<>();
        Set<String> setOfUniqueHash = new HashSet<>();
        for(SubmittedVariantEntity sve: sveList){
            String sveNewHash = mapOldHashNewHash.get(sve.getHashedMessage())
            // try to add the hash in the  set of unique hashes
            boolean res = setOfUniqueHash.add(sveNewHash)
            //if it fails to add (res=false), that means we have already seen this hash
            if (!res){
                hashCollisionSVEInBatch.add(sve);
            }
        }

        logger.info("SVE that will cause hash collision in the Batch: " + hashCollisionSVEInBatch)
        return hashCollisionSVEInBatch
    }

    void processNoHashCollision(List<SubmittedVariantEntity> noHashCollisionSVEList, Map<String,String> mapOldHashNewHash,
                                Class ssClass){
        // insert SVEs with new hash (can't update in place as _id is immutable)
        logger.info("start: insert new hash SVE : " + LocalDateTime.now())
        insertSVEWithNewHash(noHashCollisionSVEList, mapOldHashNewHash)
        logger.info("done: insert new hash SVE : " + LocalDateTime.now())

        // insert update operations for all the above inserts (SVEs with new hash) in SVOE
        logger.info("start: insert update operation into SVOE : " + LocalDateTime.now())
        insertSVOEUpdateOp(noHashCollisionSVEList, ssClass)
        logger.info("done: insert update operation into SVOE : " + LocalDateTime.now())

        // remove all impacted SVEs (SVE with old hash) from the table
        logger.info("start: remove impacted SVE : " + LocalDateTime.now())
        removeImpactedSVE(noHashCollisionSVEList, ssClass)
        logger.info("end: remove impacted SVE : " + LocalDateTime.now())
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
                                logger.info("start: processing batch : " + LocalDateTime.now())
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
                                    logger.info("start: get existing hash SVE : " + LocalDateTime.now())
                                    Set<SubmittedVariantEntity> alreadyExistingHashSVEList = getExistingHashSVEList(mapOldHashNewHash.values())
                                    Set<String> alreadyExistingHash = alreadyExistingHashSVEList.stream()
                                            .map(sve -> sve.getHashedMessage())
                                            .collect(Collectors.toSet())
                                    logger.info("done: get existing hash SVE : " + LocalDateTime.now())

                                    // partition the impacted sve list into 2 parts  (hash collision and no hash collision)
                                    Map<Boolean, List<SubmittedVariantEntity>> svePartitionMap = impactedSVEList.stream()
                                            .collect(Collectors.groupingBy { sve -> alreadyExistingHash.contains(mapOldHashNewHash.get(sve.getHashedMessage())) })

                                    // sve with no hash collision
                                    List<SubmittedVariantEntity> noHashCollisionSVEList = svePartitionMap.get(Boolean.FALSE)
                                    logger.info("Impacted sve List (No Hash Collision): " + noHashCollisionSVEList)

                                    // sve with hash collision
                                    List<SubmittedVariantEntity> hashCollisionSVEList = svePartitionMap.get(Boolean.TRUE)

                                    // check if there can be a hash collision in the batch
                                    // i.e. if there are any 2 or more SVEs that has the same new hash
                                    // when one is inserted, 2nd will fail with duplicate key error
                                    int numOfUniqueNewHash = mapOldHashNewHash.values().stream().collect(Collectors.toSet()).size()
                                    if(numOfUniqueNewHash != impactedSVEList.size()){
                                        List<SubmittedVariantEntity> hashCollisionInBatch = findHashCollisionInBatch(noHashCollisionSVEList, mapOldHashNewHash)
                                        //remove the SVEs from the no hash collision list
                                        noHashCollisionSVEList.removeAll(hashCollisionInBatch)
                                        // add the SVEs into the hash collision list
                                        hashCollisionSVEList.addAll(hashCollisionInBatch)
                                    }

                                    // process sve with no hash collision
                                    processNoHashCollision(noHashCollisionSVEList, mapOldHashNewHash, ssClass)

                                    // process sve with hash collision
                                    if (hashCollisionSVEList != null && !hashCollisionSVEList.isEmpty()) {
                                        logger.info("Impacted sve List (Hash Collision + hash collision in batch): " + hashCollisionSVEList)
                                        // capture in file - sve hash collision details
                                        logger.info("start: update hash collision file: " + LocalDateTime.now())
                                        updateFileWithHashCollisionList(hashCollisionSVEList, mapOldHashNewHash, alreadyExistingHashSVEList)
                                        logger.info("end: update hash collision file: " + LocalDateTime.now())
                                    }

                                    //update existing SVE Operations
                                    logger.info("start: update existing SVOE : " + LocalDateTime.now())
                                    updateExistingSVOE(impactedSVEList, mapOldHashNewHash)
                                    logger.info("done: update existing SVOE : " + LocalDateTime.now())
                                }

                                totalProcessed = totalProcessed + currentBatchSize
                                logger.info("total Processed till now = " + totalProcessed)
                                logger.info("total impacted till now = " + totalImpacted)

                                logger.info("done: processing batch : " + LocalDateTime.now())
                            }
                    }
                }
        }

    }
}

