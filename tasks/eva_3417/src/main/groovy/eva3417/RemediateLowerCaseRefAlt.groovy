package eva3417

import com.google.gson.Gson
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
    private Gson gson

    RemediateLowerCaseNucleotide(envPropertiesFile, hashcollisionFile) {
        this.dbEnv = createFromSpringContext(envPropertiesFile, Application.class)
        this.hashCollisionFilePath = hashcollisionFile
        this.hashingFunction = new SubmittedVariantSummaryFunction().andThen(new SHA1HashingFunction())
        this.gson = new Gson()
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

    List<SubmittedVariantEntity> getExistingHashSVEList(Collection<String> newHashList) {
        List<SubmittedVariantEntity> existingHashSVEList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(newHashList)), ssClass))
                .flatten()

        return existingHashSVEList
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
        dbEnv.mongoTemplate.insert(sveWithNewHashList, ssClass)
    }

    void processHashCollision(List<SubmittedVariantEntity> hashCollisionSVEList,
                              List<SubmittedVariantEntity> alreadyExistingHashSVEList, Class ssClass){
        // In most cases, we should have only one SVE that has a collision
        // but there can be cases where we have 2 or more SVE in the batch that has a collision
        Map<String, List<SubmittedVariantEntity>> mapOfNewHashAndSVEList = hashCollisionSVEList.stream()
                .collect(Collectors.groupingBy(sve->getSVENewHash(sve), Collectors.toList()))
        Map<String, SubmittedVariantEntity> mapOfHashAndExistingSVEInDB = alreadyExistingHashSVEList.stream()
                .collect(Collectors.toMap(sve->sve.getHashedMessage(), sve->sve))

        List<SubmittedVariantEntity> SVEToDiscardList = new ArrayList<>()
        Map<Long, List<SubmittedVariantEntity>> SVEToMergeMap =  new HashMap<>()
        List<SubmittedVariantEntity> SVEToInsertList = new ArrayList<>()
        List<SubmittedVariantEntity> SVEToRemoveList = new ArrayList<>()

        // go through each hash collision one at a time
        for(Map.Entry<String, List<SubmittedVariantEntity>> entry: mapOfNewHashAndSVEList.entrySet()){
            // all SVE with same new hash
            List<SubmittedVariantEntity> sveList = entry.getValue()
            // get the SVE in DB
            SubmittedVariantEntity sveInDB = mapOfHashAndExistingSVEInDB.get(entry.getKey())

            //assert all SVEs involved has same RS
            Set<Long> rsSet = sveList.stream().map(sve->sve.getClusteredVariantAccession())
                    .collect(Collectors.toSet())
            rsSet.add(sveInDB.getClusteredVariantAccession())

            if(rsSet.size() > 1){
                //all sve involved does not have same RS, document SVEs in file - don't process
                documentInFile(sveList)
            } else{
                //sort the SVEs involved in asc order of accession and created date
                List<SubmittedVariantEntity> sortedList = sveList.stream()
                        .sorted(Comparator.comparingLong(SubmittedVariantEntity::getAccession)
                                .thenComparing({ SubmittedVariantEntity sve -> sve.getCreatedDate() }))
                        .collect(Collectors.toList())


                // assume the one in db has the smallest accession
                Long lowestAcc = sveInDB.getAccession()
                int startIndex = 0

                SubmittedVariantEntity firstSVEInSorted = sortedList.get(0)

                // update values if the lowest accession is not in db but in sorted List
                if(firstSVEInSorted.getAccession() < sveInDB.getAccession()){
                    // we need to process first element separately (it needs to be inserted not discarded or merged)
                    startIndex = 1
                    lowestAcc = firstSVEInSorted.getAccession()
                    // Add the SVE in DB to the list for processing (will be merged)
                    sortedList.add(sveInDB)
                    // need to remove this SVE (will be inserted after fixing lower nucleotide)
                    SVEToRemoveList.add(firstSVEInSorted)
                    // need to insert the SVE from the sorted list (with lowest accession) into db
                    SVEToInsertList.add(firstSVEInSorted)
                }

                for(int i=startIndex; i<sortedList.size(); i++){
                    SubmittedVariantEntity currSVE = sortedList.get(i)
                    if(currSVE.getAccession() == lowestAcc){
                        // if acc is same - discard
                        SVEToDiscardList.add(currSVE)
                    }else{
                        // if acc is different - merge
                        SVEToMergeMap.computeIfAbsent(lowestAcc, k -> new ArrayList<>()).add(currSVE)
                    }
                }
            }
        }

        // discard SVEs
        logger.info("Hash Collision Discard SVE List: " + gson.toJson(SVEToDiscardList))
        discardSVEAndInsertOperations(SVEToDiscardList, ssClass)

        // merge SVEs
        mergeSVEAndInsertOperations(SVEToMergeMap, ssClass)

        // Insert SVEs
        logger.info("Hash Collision Insert SVE List: " + gson.toJson(SVEToInsertList))
        insertSVEWithNewHash(SVEToInsertList, ssClass)

        // delete the corrected SVEs
        logger.info("Hash Collision Remove SVE List: " + gson.toJson(SVEToRemoveList))
        removeImpactedSVE(SVEToRemoveList, ssClass)

    }

    void discardSVEAndInsertOperations(List<SubmittedVariantEntity> discardSVEList, Class ssClass){
        List<SubmittedVariantOperationEntity> svoeList = new ArrayList<>()
        for(SubmittedVariantEntity sve: discardSVEList){
            // we need to still correct the lower case issue before we put this sve in inactiveObjects of the operation
            SubmittedVariant sveModel = sve.getModel()
            sveModel.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
            sveModel.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
            SubmittedVariantEntity submittedVariantEntity = new SubmittedVariantEntity(sve.getAccession(),
                    getSVENewHash(sve), sveModel, sve.getVersion())

            // create discard operation
            SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity()
            svoe.fill(EventType.DISCARDED, sve.getAccession(), "Submitted variant discarded due to duplicate hash",
                    Collections.singletonList(new SubmittedVariantInactiveEntity(submittedVariantEntity)))
            svoe.setId("EVA3417_DISCARD_SS_"+sve.getAccession()+"_HASH_" + sve.getHashedMessage())

            //add to the list of operations
            svoeList.add(svoe)
        }

        // delete the discarded SVEs
        removeImpactedSVE(discardSVEList, ssClass)

        // Add the discard operation for all SVEs
        dbEnv.mongoTemplate.insert(svoeList, ssClass.equals(sveClass) ? svoeClass : dbsnpSvoeClass)
    }


    void mergeSVEAndInsertOperations(Map<Long, List<SubmittedVariantEntity>> mergeSVEMap, Class ssClass){
        List<SubmittedVariantOperationEntity> svoeList = new ArrayList<>()
        for(Map.Entry<Long, List<SubmittedVariantEntity>> entry: mergeSVEMap){
            Long lowestAcc = entry.getKey()
            for(SubmittedVariantEntity sve: entry.getValue()) {
                // we need to still correct the lower case issue before we put this sve in inactiveObjects of the operation
                SubmittedVariant sveModel = sve.getModel()
                sveModel.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
                sveModel.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
                SubmittedVariantEntity submittedVariantEntity = new SubmittedVariantEntity(sve.getAccession(),
                        getSVENewHash(sve), sveModel, sve.getVersion())

                // create merge operation
                SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity()
                svoe.fill(EventType.MERGED, sve.getAccession(), lowestAcc,
                        "After fixing lowercase nucleotide issue, variant merged due to duplicate hash",
                        Collections.singletonList(new SubmittedVariantInactiveEntity(submittedVariantEntity)))
                svoe.setId("EVA3417_MERGED_"+sve.getAccession()+"_INTO_"+lowestAcc+"_HASH_" + sve.getHashedMessage() )

                //add to the list of operations
                svoeList.add(svoe)
            }
        }

        // delete the merged SVEs
        List<SubmittedVariantEntity> mergeSVEList = mergeSVEMap.entrySet().stream()
                .flatMap(entry->entry.getValue().stream())
                .collect(Collectors.toList())
        logger.info("Hash Collision Merge SVE List: " + gson.toJson(mergeSVEList))
        removeMergedSVE(mergeSVEList, ssClass)

        // Add the merge operation for all SVEs
        dbEnv.mongoTemplate.insert(svoeList, ssClass.equals(sveClass) ? svoeClass : dbsnpSvoeClass)
    }

    void removeMergedSVE(List<SubmittedVariantEntity> mergeSVEList, Class ssClass){
        List<String> mergedSVEIds = mergeSVEList.stream()
                .map(sve -> sve.getHashedMessage())
                .collect(Collectors.toList())
        List<String> mergedSVEAcc = mergeSVEList.stream()
                .map(sve -> sve.getAccession())
                .collect(Collectors.toList())

        Query findQuery = query(where("_id").in(mergedSVEIds).and("accession").in(mergedSVEAcc))

        dbEnv.mongoTemplate.remove(findQuery, ssClass)
    }


    void documentInFile(List<SubmittedVariantEntity> sveList) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(this.hashCollisionFilePath, true))) {
            writer.write(sveList)
        } catch (IOException e) {
            e.printStackTrace()
        }
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
                .and("_id").not().regex("^EVA3417_UPDATED_")
        new RetryableBatchingCursor<SubmittedVariantOperationEntity>(filterCriteria, dbEnv.mongoTemplate, ssClass).each {
            sveOperationBatchList ->
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


    void insertSVOEUpdateOp(List<SubmittedVariantEntity> sveList, Class ssClass) {
        def sveUpdateOpList = new ArrayList<>()
        for (SubmittedVariantEntity sve : sveList) {
            def updateOpToWrite = ssClass.equals(sveClass)
                    ? new SubmittedVariantOperationEntity() : new DbsnpSubmittedVariantOperationEntity()
            updateOpToWrite.fill(EventType.UPDATED, sve.getAccession(), null,
                    "EVA3417 - SS updated with upper case ref, alt and a new hash",
                    Arrays<SubmittedVariantInactiveEntity>.asList(new SubmittedVariantInactiveEntity(sve)).asList())

            updateOpToWrite.setId("EVA3417_UPDATED_"+sve.getReferenceSequenceAccession()+"_"+sve.getAccession())

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


    void processNoHashCollision(List<SubmittedVariantEntity> noHashCollisionSVEList, Class ssClass){
        // insert SVEs with new hash (can't update in place as _id is immutable)
        insertSVEWithNewHash(noHashCollisionSVEList, ssClass)

        // insert update operations for all the above inserts (SVEs with new hash) in SVOE
        insertSVOEUpdateOp(noHashCollisionSVEList, ssClass)

        // remove all impacted SVEs (SVE with old hash) from the table
        removeImpactedSVE(noHashCollisionSVEList, ssClass)
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
                                    List<SubmittedVariantEntity> alreadyExistingSVEList = getExistingHashSVEList(mapOldHashNewHash.values())
                                    Set<String> alreadyExistingHash = alreadyExistingSVEList.stream()
                                            .map(sve -> sve.getHashedMessage())
                                            .collect(Collectors.toSet())

                                    // partition the impacted sve list into 2 parts (hash collision and no hash collision)
                                    Map<Boolean, List<SubmittedVariantEntity>> svePartitionMap = impactedSVEList.stream()
                                            .collect(Collectors.groupingBy(sve -> alreadyExistingHash.contains(mapOldHashNewHash.get(sve.getHashedMessage()))))

                                    // sve with no hash collision
                                    List<SubmittedVariantEntity> noHashCollisionSVEList = svePartitionMap.get(Boolean.FALSE)
                                    if(noHashCollisionSVEList == null){
                                        noHashCollisionSVEList = new ArrayList<>()
                                    }

                                    // sve with hash collision
                                    List<SubmittedVariantEntity> hashCollisionSVEList = svePartitionMap.get(Boolean.TRUE)
                                    if(hashCollisionSVEList == null){
                                        hashCollisionSVEList = new ArrayList<>()
                                    }

                                    // check if there can be a hash collision in the batch
                                    // i.e. if there are any 2 or more SVEs that has the same new hash
                                    int numOfUniqueNewHash = mapOldHashNewHash.values().stream().collect(Collectors.toSet()).size()
                                    if(numOfUniqueNewHash != impactedSVEList.size()){
                                        // go through the no hash collision list and make a map of new hash and list of SVE with that new hash
                                        // then filter and take only those where the list has more than one SVE
                                        List<Map.Entry<String, List<SubmittedVariantEntity>>> mapEntryOfNewHashAndSVEList = noHashCollisionSVEList.stream()
                                                .collect(Collectors.groupingBy(sve->getSVENewHash(sve), Collectors.toList()))
                                                .entrySet().stream().filter(entry->entry.getValue().size() > 1)
                                        .collect(Collectors.toList())

                                        // for each hash, keep the first sve in the no hash collision list and move the others to hash collision list
                                        for(Map.Entry<String, List<SubmittedVariantEntity>> entry: mapEntryOfNewHashAndSVEList){
                                            List<SubmittedVariantEntity> sveList = entry.getValue()
                                            for(int i=1;i<sveList.size();i++){
                                                // remove from the no hash collision list
                                                noHashCollisionSVEList.remove(sveList.get(i))
                                                // add to the hash collision list
                                                hashCollisionSVEList.add(sveList.get(i))
                                            }

                                            // add the first SVE in the list to the alreadyExisting SVE list,
                                            // an equivalent of this will be inserted in db as part of no hash collision list
                                            // before we start processing hash collision list
                                            SubmittedVariantEntity firstSVE = sveList.get(0)
                                            SubmittedVariant sveModel = firstSVE.getModel()
                                            sveModel.setReferenceAllele(firstSVE.getReferenceAllele().toUpperCase())
                                            sveModel.setAlternateAllele(firstSVE.getAlternateAllele().toUpperCase())
                                            SubmittedVariantEntity submittedVariantEntity = new SubmittedVariantEntity(firstSVE.getAccession(),
                                                    getSVENewHash(firstSVE), sveModel, firstSVE.getVersion())

                                            alreadyExistingSVEList.add(submittedVariantEntity)
                                        }
                                    }

                                    // process sve with no hash collision
                                    logger.info("Impacted sve List (No Hash Collision): " + gson.toJson(noHashCollisionSVEList))
                                    processNoHashCollision(noHashCollisionSVEList, ssClass)

                                    // process sve with hash collision
                                    if (hashCollisionSVEList != null && !hashCollisionSVEList.isEmpty()) {
                                        logger.info("Impacted sve List (Hash Collision + hash collision in batch): " + gson.toJson(hashCollisionSVEList))
                                        processHashCollision(hashCollisionSVEList, alreadyExistingSVEList, ssClass)
                                    }

                                    // update existing SVE Operations
                                    logger.info("start: update existing SVOE : " + LocalDateTime.now())
                                    updateExistingSVOE(impactedSVEList, mapOldHashNewHash, ssClass.equals(sveClass) ? svoeClass : dbsnpSvoeClass)
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

