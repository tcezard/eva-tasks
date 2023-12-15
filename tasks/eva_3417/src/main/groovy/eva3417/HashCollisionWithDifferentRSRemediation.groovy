package eva3417

import com.google.gson.Gson
import com.google.gson.reflect.TypeToken
import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.BulkOperations
import org.springframework.data.mongodb.core.query.Query
import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.summary.SubmittedVariantSummaryFunction
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment

import java.nio.file.Files
import java.util.function.Function
import java.util.stream.Collectors

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

def cli = new CliBuilder()
cli.envPropertiesFile(args: 1, "properties file for connecting to db", required: true)
cli.baseDirPath(args: 2, "Path to the base dir of logs", required: true)
cli.summaryFilePath(args: 3, "Path to the summary file dir", required: true)

def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

new SummariseHashCollision(options.envPropertiesFile, options.baseDirPath, options.summaryFilePath).runSummary()

class SummariseHashCollision {
    private EVADatabaseEnvironment dbEnv
    private String baseDirPath
    private String summaryFilePath
    private Function<ISubmittedVariant, String> hashingFunction
    static def logger = LoggerFactory.getLogger(Application.class)
    private Gson gson

    private Map<String, List<CollisionSummary>> collisionSummaryMap = new HashMap<>()

    SummariseHashCollision(envPropertiesFile, baseDirPath, summaryFilePath) {
        this.dbEnv = createFromSpringContext(envPropertiesFile, Application.class)
        this.baseDirPath = baseDirPath
        this.summaryFilePath = summaryFilePath
        this.hashingFunction = new SubmittedVariantSummaryFunction().andThen(new SHA1HashingFunction())
        this.gson = new Gson()
    }

    String getSVENewHash(SubmittedVariantEntity sve) {
        SubmittedVariant temp = sve.getModel()
        temp.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
        temp.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
        return hashingFunction.apply(temp)
    }

    ArrayList<SubmittedVariantEntity> getListOfSVEFromLine(String sveListString) {
        ArrayList<SubmittedVariantEntity> sveList = (ArrayList<SubmittedVariantEntity>) gson.fromJson(sveListString,
                new TypeToken<ArrayList<SubmittedVariantEntity>>() {}.getType())

        return sveList
    }

    Map<String, SubmittedVariantEntity> getAllSVEInDBForNewHashes(List<List<SubmittedVariantEntity>> sveInCollisionFile) {
        // get all new hash sves
        Set<String> setOfHashes = sveInCollisionFile.stream()
                .flatMap(sveList -> sveList.stream())
                .map(sve -> getSVENewHash(sve))
                .collect(Collectors.toSet())

        List<SubmittedVariantEntity> dbSVEList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(setOfHashes)), ssClass))
                .flatten()

        return dbSVEList.stream().collect(Collectors.toMap(sve -> sve.getHashedMessage(), sve -> sve))
    }

    Map<Long, List<ClusteredVariantEntity>> getAllCVEInDB(List<List<SubmittedVariantEntity>> sveInCollisionFile,
                                                          Map<String, SubmittedVariantEntity> sveInDBMap) {
        Set<ClusteredVariantEntity> setOfAccessions = new HashSet<>()
        // get accesison from sves in file
        Set<Long> accessionsFromFileSVE = sveInCollisionFile.stream()
                .flatMap(sveList -> sveList.stream())
                .filter(sve -> sve.getClusteredVariantAccession() != null)
                .map(sve -> sve.getClusteredVariantAccession())
                .collect(Collectors.toSet())
        setOfAccessions.addAll(accessionsFromFileSVE)

        // get accessions from db
        Set<Long> accessionsFromDBSVE = sveInDBMap.values().stream()
                .filter(sve -> sve.getClusteredVariantAccession() != null)
                .map(sve -> sve.getClusteredVariantAccession())
                .collect(Collectors.toSet())
        // add accessions from both the sve in file as well as db
        setOfAccessions.addAll(accessionsFromDBSVE)

        List<ClusteredVariantEntity> cveInDB = [cveClass, dbsnpCveClass]
                .collectMany(cveClass -> dbEnv.mongoTemplate.find(query(where("accession").in(setOfAccessions)), cveClass))
                .flatten()

        Map<Long, List<ClusteredVariantEntity>> cveInDBMap = cveInDB.stream().collect(Collectors.groupingBy(cve -> cve.getAccession()))

        return cveInDBMap
    }

    void buildSummary() {
        String[] affectedAsmList = [
                "GCA_000002285.2", "GCA_000002315.3", "GCA_000003055.5", "GCA_011100555.1", "GCA_000002035.3",
                "GCA_000001635.6", "GCA_000349105.1", "GCA_016772045.1", "GCA_000003025.6", "GCA_000349185.1",
                "GCA_000247795.2", "GCA_000473445.2",
                "GCA_001433935.1"
        ]

        for (String asm : affectedAsmList) {
            File asmHashCollisionFile = new File(baseDirPath + "/" + asm + "/hash_collision.txt")
            // read all the hash collision sve from the list
            List<List<SubmittedVariantEntity>> sveInCollisionFile = Files.lines(asmHashCollisionFile.toPath())
                    .map(line -> getListOfSVEFromLine(line))
                    .collect(Collectors.toList())

            // Get all SVE in DB For New Hashes
            Map<String, SubmittedVariantEntity> sveInDBMap = getAllSVEInDBForNewHashes(sveInCollisionFile)
            // Get all RS in DB for SVE in file and SVE in DB
            Map<Long, List<ClusteredVariantEntity>> cveInDB = getAllCVEInDB(sveInCollisionFile, sveInDBMap)

            collisionSummaryMap.put(asm, new ArrayList<CollisionSummary>())

            for (List<SubmittedVariantEntity> sveList : sveInCollisionFile) {
                if (sveList.size() != 1) {
                    // TODO: only one variant that has this condition, needs to be checked
                    logger.error("skip processing as svelist size is not equal to 1 ")
                    continue
                }
                SubmittedVariantEntity sveInFile = sveList.get(0)
                SubmittedVariantEntity sveInDB = sveInDBMap.get(getSVENewHash(sveInFile))
                collisionSummaryMap.get(asm).add(new CollisionSummary(sveInFile, sveInDB, cveInDB))
            }
        }
    }

    void writeSummaryToFile() {
        try (BufferedWriter br = new BufferedWriter(new FileWriter(summaryFilePath))) {
            String[] header = new String[]{"assembly",
                    "sve_in_file", "file_accession", "file_rs", "file_rs_valid", "file_remapped_from",
                    "sve_in_db", "db_accession", "db_rs", "db_rs_valid", "db_remapped_from",
                    "file_ref", "file_alt", "db_ref", "db_alt",
                    "category"}
            br.write(header.join(","))
            br.write("\n")

            for (Map.Entry<String, List<CollisionSummary>> asmSummary : collisionSummaryMap.entrySet()) {
                for (CollisionSummary summary : asmSummary.getValue()) {
                    String[] row = new String[]{
                            asmSummary.getKey().toString(),

                            summary.getSveInFile().getHashedMessage(), summary.getSveInFile().getAccession(),
                            summary.getSveInFile().getClusteredVariantAccession(), summary.getRsInFileValid(),
                            summary.getSveInFile().getRemappedFrom(),

                            summary.getSveInDB().getHashedMessage(), summary.getSveInDB().getAccession(),
                            summary.getSveInDB().getClusteredVariantAccession(), summary.getRsInDBValid(),
                            summary.getSveInDB().getRemappedFrom(),

                            summary.getSveInFile().getReferenceAllele(), summary.getSveInFile().getAlternateAllele(),
                            summary.getSveInDB().getReferenceAllele(), summary.getSveInDB().getAlternateAllele(),

                            summary.getCategory()}
                    br.write(row.join(","))
                    br.write("\n")
                }
            }
        }
    }

    void remediateCategoryOneValidRS() {
        logger.info("Processing SVEs where only one of them has a valid RS ")

        for (Map.Entry<String, List<CollisionSummary>> asmSummary : collisionSummaryMap.entrySet()) {
            String assembly = asmSummary.getKey()
            logger.info("processing for assembly : " + assembly)
            List<CollisionSummary> collisionSummaryList = asmSummary.getValue()
            // filter and take sve where we have  at the most one valid RS (in this method we will remediate only these)
            List<CollisionSummary> listOfValidCollision = collisionSummaryList.stream()
                    .filter(colsum -> colsum.getCategory() != CATEGORY.BOTH_SVE_HAS_RS_BOTH_VALID)
                    .collect(Collectors.toList())
            logger.info("Number of collisions with no or one valid rs: " + listOfValidCollision.size())
            if (listOfValidCollision == null || listOfValidCollision.isEmpty()) {
                logger.info("No collision with 0 or 1 valid RS found in assembly " + assembly)
                continue
            }

            // find out the residing collection for both sve in file and sve in db
            Set<String> sveInFileHashes = listOfValidCollision.stream()
                    .map(colsum -> colsum.getSveInFile().getHashedMessage())
                    .collect(Collectors.toSet())
            Set<String> sveInDBHashes = listOfValidCollision.stream()
                    .map(colsum -> colsum.getSveInDB().getHashedMessage())
                    .collect(Collectors.toSet())
            Set<String> allHashes = new HashSet<>()
            allHashes.addAll(sveInFileHashes)
            allHashes.addAll(sveInDBHashes)

            Set<String> sveCollection = dbEnv.mongoTemplate.find(query(where("_id").in(allHashes)), sveClass)
                    .stream().map(sve -> sve.getHashedMessage())
                    .collect(Collectors.toSet())
            Set<String> dbsnpSVECollection = dbEnv.mongoTemplate.find(query(where("_id").in(allHashes)), dbsnpSveClass)
                    .stream().map(sve -> sve.getHashedMessage())
                    .collect(Collectors.toSet())

            List<SubmittedVariantEntity> sveDeleteList = new ArrayList<>()
            List<SubmittedVariantEntity> sveInsertList = new ArrayList<>()
            List<SubmittedVariantOperationEntity> svoeInsertList = new ArrayList<>()
            List<SubmittedVariantEntity> sveRSUpdateList = new ArrayList<>()
            List<SubmittedVariantEntity> dbsnpSveDeleteList = new ArrayList<>()
            List<SubmittedVariantEntity> dbsnpSveInsertList = new ArrayList<>()
            List<SubmittedVariantOperationEntity> dbsnpSvoeInsertList = new ArrayList<>()
            List<SubmittedVariantEntity> dbsnpSveRSUpdateList = new ArrayList<>()

            for (CollisionSummary colSummary : listOfValidCollision) {
                //  There can be different cases here
                //  1. if both sve are not remapped
                //      - take the one with lowest created date and accession (copy rs from other sve if required)
                //  2. if both sve are remapped
                //      - see what sve accession exist in the original assembly
                //          - if only one exist, take that one
                //          - if both exist or both do not exist, take the one with lowest created date/accession
                //  3. if only one is remapped
                //      - take the non-remapped one
                // in all cases, if the one we are keeping does not have a valid rs, copy valid rs from the other

                SubmittedVariantEntity sveInFile = colSummary.getSveInFile()
                SubmittedVariantEntity sveInDB = colSummary.getSveInDB()
                // assume we need to keep the sve in db and merge the sve in file (lowercase)
                SubmittedVariantEntity sveToKeep = sveInDB
                boolean sveToKeepRSValid = colSummary.getRsInDBValid()
                SubmittedVariantEntity sveToMerge = sveInFile
                boolean sveToMergeRSValid = colSummary.getRsInFileValid()
                Boolean sveToKeepIsInFile = false

                if (sveInFile.getRemappedFrom() == null && sveInDB.getRemappedFrom() == null) {
                    // case 1: both are not remapped (take the one with lowest created date/accession)
                    if (sveInFile.getCreatedDate().isBefore(sveInDB.getCreatedDate())) {
                        sveToKeep = sveInFile
                        sveToMerge = sveInDB
                        sveToKeepRSValid = colSummary.getRsInFileValid()
                        sveToMergeRSValid = colSummary.getRsInDBValid()
                        sveToKeepIsInFile = true
                    } else if (sveInFile.getCreatedDate().isEqual(sveInDB.getCreatedDate())) {
                        if (sveInFile.getAccession() < sveInDB.getAccession()) {
                            sveToKeep = sveInFile
                            sveToMerge = sveInDB
                            sveToKeepRSValid = colSummary.getRsInFileValid()
                            sveToMergeRSValid = colSummary.getRsInDBValid()
                            sveToKeepIsInFile = true
                        }
                    }
                } else if (sveInFile.getRemappedFrom() != null && sveInDB.getRemappedFrom() != null) {
                    // case 2: both are remapped check which sve accession is present in original assembly that gets the priority
                    List<SubmittedVariantEntity> sveInFileOrgAsm = findSVEAccessionInOrgAssembly(sveInFile.getRemappedFrom(), sveInFile.getAccession())
                    List<SubmittedVariantEntity> sveInDBOrgAsm = findSVEAccessionInOrgAssembly(sveInDB.getRemappedFrom(), sveInDB.getAccession())
                    if (((sveInFileOrgAsm != null && !sveInFileOrgAsm.isEmpty()) && (sveInDBOrgAsm != null && !sveInDBOrgAsm.isEmpty()))
                            || ((sveInFileOrgAsm == null || sveInFileOrgAsm.isEmpty()) && (sveInDBOrgAsm == null || sveInDBOrgAsm.isEmpty()))) {
                        // both sve accession are present in original assembly or both sve accession are not present in the original assembly
                        // similar to case 1 - take the one with lower created date and sve accession
                        if (sveInFile.getCreatedDate().isBefore(sveInDB.getCreatedDate())) {
                            sveToKeep = sveInFile
                            sveToMerge = sveInDB
                            sveToKeepRSValid = colSummary.getRsInFileValid()
                            sveToMergeRSValid = colSummary.getRsInDBValid()
                            sveToKeepIsInFile = true
                        } else if (sveInFile.getCreatedDate().isEqual(sveInDB.getCreatedDate())) {
                            if (sveInFile.getAccession() < sveInDB.getAccession()) {
                                sveToKeep = sveInFile
                                sveToMerge = sveInDB
                                sveToKeepRSValid = colSummary.getRsInFileValid()
                                sveToMergeRSValid = colSummary.getRsInDBValid()
                                sveToKeepIsInFile = true
                            }
                        }
                    } else {
                        // only one asm accession is present in original assembly, take the one that is present
                        if (sveInDBOrgAsm == null || sveInDBOrgAsm.isEmpty()) {
                            sveToKeep = sveInFile
                            sveToMerge = sveInDB
                            sveToKeepRSValid = colSummary.getRsInFileValid()
                            sveToMergeRSValid = colSummary.getRsInDBValid()
                            sveToKeepIsInFile = true
                        }
                    }
                } else if (sveInFile.getRemappedFrom() == null || sveInDB.getRemappedFrom() == null) {
                    // case 3: only one is remapped , non remapped one gets the priority
                    if (sveInFile.getRemappedFrom() == null) {
                        sveToKeep = sveInFile
                        sveToMerge = sveInDB
                        sveToKeepRSValid = colSummary.getRsInFileValid()
                        sveToMergeRSValid = colSummary.getRsInDBValid()
                        sveToKeepIsInFile = true
                    }
                }

                Class sveToKeepCollection = getCollection(sveCollection, dbsnpSVECollection, sveToKeep)
                Class sveToMergeCollection = getCollection(sveCollection, dbsnpSVECollection, sveToMerge)

                // create a merge operation for sve to merge and delete the merged sve
                if (sveToMergeCollection == sveClass) {
                    svoeInsertList.add(getSVOEForMergeOperation(sveToMerge, sveToKeep.getAccession()))
                    sveDeleteList.add(sveToMerge)
                } else {
                    dbsnpSvoeInsertList.add(getSVOEForMergeOperation(sveToMerge, sveToKeep.getAccession()))
                    dbsnpSveDeleteList.add(sveToMerge)
                }

                // if sve to keep does not have a valid rs but sve to merge has, copy it to sve to keep
                if (sveToKeepRSValid==false && sveToMergeRSValid==true) {
                    sveToKeep.setClusteredVariantAccession(sveToMerge.getClusteredVariantAccession())
                }

                if (sveToKeepIsInFile) {
                    // delete the lowercase sve (will insert the remediated one)
                    sveDeleteList.add(sveToKeep)
                    // remediate and add the sve to the db
                    if (sveToKeepCollection == sveClass) {
                        sveInsertList.add(getRemediatedSVE(sveToKeep))
                    } else {
                        dbsnpSveInsertList.add(getRemediatedSVE(sveToKeep))
                    }
                } else {
                    // if sve to keep is present in db but does not have a valid rs, then update the rs
                    if (sveToKeepRSValid==false && sveToMergeRSValid==true) {
                        if (sveToKeepCollection == sveClass) {
                            sveRSUpdateList.add(sveToKeep)
                        } else {
                            dbsnpSveRSUpdateList.add(sveToKeep)
                        }
                    }
                }
            }

            // insert the merge operation entries for all SVEs
            logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Merged Operations List (SVOE): " + gson.toJson(svoeInsertList))
            logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Merged Operations List (DBSNPSVOE): " + gson.toJson(dbsnpSvoeInsertList))
            if (svoeInsertList != null && !svoeInsertList.isEmpty()) {
                dbEnv.mongoTemplate.insert(svoeInsertList, svoeClass)
            }
            if (dbsnpSvoeInsertList != null && !dbsnpSvoeInsertList.isEmpty()) {
                dbEnv.mongoTemplate.insert(dbsnpSvoeInsertList, dbsnpSvoeClass)
            }

            // delete the merged SVEs - delete first before inserting in order to avoid deleting the remediated one
            logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Delete Merged SVE : " + gson.toJson(sveDeleteList))
            if (sveDeleteList != null && !sveDeleteList.isEmpty()) {
                removeMergedSVE(sveDeleteList, sveClass)
            }
            if (dbsnpSveDeleteList != null && !dbsnpSveDeleteList.isEmpty()) {
                removeMergedSVE(dbsnpSveDeleteList, dbsnpSveClass)
            }

            // insert updated SVEs into DB
            logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Insert SVE : " + gson.toJson(sveInsertList))
            logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Insert DBSNPSVE: " + gson.toJson(dbsnpSveInsertList))
            if (sveInsertList != null && !sveInsertList.isEmpty()) {
                dbEnv.mongoTemplate.insert(sveInsertList, sveClass)
            }
            if (dbsnpSveInsertList != null && !dbsnpSveInsertList.isEmpty()) {
                dbEnv.mongoTemplate.insert(dbsnpSveInsertList, dbsnpSveClass)
            }

            if (!sveRSUpdateList.isEmpty() || !dbsnpSveRSUpdateList.isEmpty()) {
                updateRSForSVe(sveRSUpdateList, dbsnpSveRSUpdateList)
            }
        }
    }

    void updateRSForSVe(List<SubmittedVariantEntity> sveUpdateList, List<SubmittedVariantEntity> dbsnpSveUpdateList) {
        logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Update RS in SVE : " + gson.toJson(sveUpdateList))
        logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Update RS in DBSNPSVE : " + gson.toJson(dbsnpSveUpdateList))

        def sveBulkUpdates = dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, sveClass)
        def dbsnpSveBulkUpdates = dbEnv.mongoTemplate.bulkOps(BulkOperations.BulkMode.UNORDERED, dbsnpSveClass)

        for (SubmittedVariantEntity sve : sveUpdateList) {
            Query findQuery = query(where("_id").is(sve.getHashedMessage()).and("rs").exists(false))
            Update updateQuery = new Update().set("rs", sve.getClusteredVariantAccession())
            sveBulkUpdates.updateOne(findQuery, updateQuery)
        }

        for (SubmittedVariantEntity dbsnpSve : dbsnpSveUpdateList) {
            Query findQuery = query(where("_id").is(dbsnpSve.getHashedMessage()).and("rs").exists(false))
            Update updateQuery = new Update().set("rs", dbsnpSve.getClusteredVariantAccession())
            dbsnpSveBulkUpdates.updateOne(findQuery, updateQuery)
        }

        //execute all bulk updates
        sveBulkUpdates.execute()
        dbsnpSveBulkUpdates.execute()
    }

    Class getCollection(Set<String> sveCollection, Set<String> dbsnpSVECollection, SubmittedVariantEntity sve) {
        if (sveCollection.contains(sve.getHashedMessage())) {
            return sveClass
        } else if (dbsnpSVECollection.contains(sve.getHashedMessage())) {
            return dbsnpSveClass
        } else {
            throw new RuntimeException("Could not find the collection for sve " + gson.toJson(sve))
        }
    }

    List<SubmittedVariantEntity> findSVEAccessionInOrgAssembly(String asm, Long acc) {
        List<SubmittedVariantEntity> dbSVEList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("seq").is(asm).and("accession").is(acc)), ssClass))
                .flatten()
        return dbSVEList
    }

    void removeMergedSVE(List<SubmittedVariantEntity> mergeSVEList, Class ssClass) {
        List<String> mergedSVEIds = mergeSVEList.stream()
                .map(sve -> sve.getHashedMessage())
                .collect(Collectors.toList())
        List<String> mergedSVEAcc = mergeSVEList.stream()
                .map(sve -> sve.getAccession())
                .collect(Collectors.toList())

        Query findQuery = query(where("_id").in(mergedSVEIds).and("accession").in(mergedSVEAcc))

        dbEnv.mongoTemplate.remove(findQuery, ssClass)
    }

    SubmittedVariantEntity getRemediatedSVE(SubmittedVariantEntity submittedVariantEntity) {
        SubmittedVariant sveModel = submittedVariantEntity.getModel()
        sveModel.setReferenceAllele(submittedVariantEntity.getReferenceAllele().toUpperCase())
        sveModel.setAlternateAllele(submittedVariantEntity.getAlternateAllele().toUpperCase())
        return new SubmittedVariantEntity(submittedVariantEntity.getAccession(),
                getSVENewHash(submittedVariantEntity), sveModel, submittedVariantEntity.getVersion())
    }

    SubmittedVariantOperationEntity getSVOEForMergeOperation(SubmittedVariantEntity submittedVariantEntity, Long accessionMergeInto) {
        // remediate sve if it contains lowercase nucleotide before putting in svoe
        SubmittedVariantEntity sve = getRemediatedSVE(submittedVariantEntity)

        // create merge operation
        SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity()
        svoe.fill(EventType.MERGED, sve.getAccession(), accessionMergeInto,
                "After fixing lowercase nucleotide issue, variant merged due to duplicate hash",
                Collections.singletonList(new SubmittedVariantInactiveEntity(sve)))
        svoe.setId("EVA3417_MERGED_" + sve.getAccession() + "_INTO_" + accessionMergeInto + "_HASH_" + sve.getHashedMessage())

        return svoe
    }

    void runSummary() {
        buildSummary()
        //printSummary()
        writeSummaryToFile()

        // remediate the case where only one sve has a valid rs
        remediateCategoryOneValidRS()
    }
}


enum CATEGORY {
    NO_SVE_HAS_RS,              // 0
    ONE_SVE_HAS_RS_NONE_VALID,  // 0
    ONE_SVE_HAS_RS_ONE_VALID,   // 23262
    BOTH_SVE_HAS_RS_NONE_VALID, // 51
    BOTH_SVE_HAS_RS_ONE_VALID,  // 26
    BOTH_SVE_HAS_RS_BOTH_VALID  // 161
}

class CollisionSummary {
    private SubmittedVariantEntity sveInFile
    private Long rsInFile
    private boolean rsInFileValid

    private SubmittedVariantEntity sveInDB
    private Long rsInDB
    private boolean rsInDBValid

    private CATEGORY category

    CollisionSummary(SubmittedVariantEntity sveInFile, SubmittedVariantEntity sveInDB, Map<Long, List<ClusteredVariantEntity>> cveInDB) {
        this.sveInFile = sveInFile
        this.rsInFile = sveInFile.getClusteredVariantAccession()
        this.rsInFileValid = checkIfValidRS(cveInDB, rsInFile, sveInFile.getReferenceSequenceAccession())

        this.sveInDB = sveInDB
        this.rsInDB = sveInDB.getClusteredVariantAccession()
        this.rsInDBValid = checkIfValidRS(cveInDB, rsInDB, sveInDB.getReferenceSequenceAccession())

        this.category = identifyCategory()
    }

    boolean checkIfValidRS(Map<Long, List<ClusteredVariantEntity>> cveInDB, Long rs, String asm) {
        if (rs == null) {
            return false
        } else {
            List<ClusteredVariantEntity> rsList = cveInDB.get(rs)
            if (rsList == null || rsList.isEmpty()) {
                return false
            } else {
                Optional<ClusteredVariantEntity> opRS = rsList.stream()
                        .filter(cve -> cve.getAccession() == rs && cve.getAssemblyAccession() == asm)
                        .findAny()
                if (opRS.isPresent()) {
                    return true
                } else {
                    return false
                }
            }
        }
    }

    CATEGORY identifyCategory() {
        if (rsInFile == null && rsInDB == null) {
            // none of the SVE contains any RS
            return CATEGORY.NO_SVE_HAS_RS
        } else if (rsInFile != null && rsInDB != null) {
            // both the SVE contains RS
            if (!rsInFileValid && !rsInDBValid) {
                // both the RS are invalid
                return CATEGORY.BOTH_SVE_HAS_RS_NONE_VALID
            } else if (rsInFileValid && rsInDBValid) {
                //both the RS are valid
                return CATEGORY.BOTH_SVE_HAS_RS_BOTH_VALID
            } else if (!rsInFileValid || !rsInDBValid) {
                // one RS is valid while the other is Invalid
                return CATEGORY.BOTH_SVE_HAS_RS_ONE_VALID
            }
        } else if (rsInFile != null || rsInDB != null) {
            // one of the SVE contains RS
            if (!rsInFileValid && !rsInDBValid) {
                // one of the SVE contains RS but that RS is not valid
                return CATEGORY.ONE_SVE_HAS_RS_NONE_VALID
            } else if (rsInFileValid || rsInDBValid) {
                // one of the SVE contains a valid RS
                return CATEGORY.ONE_SVE_HAS_RS_ONE_VALID
            }
        }
    }

    SubmittedVariantEntity getSveInFile() {
        return sveInFile
    }

    boolean getRsInFileValid() {
        return rsInFileValid
    }

    SubmittedVariantEntity getSveInDB() {
        return sveInDB
    }

    boolean getRsInDBValid() {
        return rsInDBValid
    }

    CATEGORY getCategory() {
        return category
    }
}


