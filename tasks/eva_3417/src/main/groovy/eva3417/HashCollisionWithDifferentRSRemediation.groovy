package eva3417

import com.google.gson.Gson
import com.google.gson.reflect.TypeToken
import groovy.cli.picocli.CliBuilder
import org.apache.commons.lang3.tuple.ImmutablePair
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Query
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
cli.summaryFilePath(args: 2, "Path to the base dir of logs", required: true)

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
        Set<ClusteredVariantEntity> accessionsFromFileSVE = sveInCollisionFile.stream()
                .flatMap(sveList -> sveList.stream())
                .filter(sve -> sve.getClusteredVariantAccession() != null)
                .map(sve -> sve.getClusteredVariantAccession())
                .collect(Collectors.toSet())
        setOfAccessions.addAll(accessionsFromFileSVE)

        // get accessions from db
        Set<ClusteredVariantEntity> accessionsFromDBSVE = sveInDBMap.values().stream()
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
                "GCA_000247795.2", "GCA_000473445.2", "GCA_001433935.1"
        ]

        for (String asm : affectedAsmList) {
            File asmHashCollisionFile = new File(baseDirPath + "/" + asm + "/hash_collision.txt")
            // read all the hash collision sve from the list (filter out if already remediated)
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
                    "sve_in_file", "file_ref", "file_alt", "file_rs", "file_rs_valid", "file_remapped_from",
                    "sve_in_db", "db_ref", "db_alt", "db_rs", "db_rs_valid", "db_remapped_from",
                    "category"}
            br.write(header.join(","))
            br.write("\n")

            for (Map.Entry<String, List<CollisionSummary>> asmSummary : collisionSummaryMap.entrySet()) {
                for (CollisionSummary summary : asmSummary.getValue()) {
                    String[] row = new String[]{
                            asmSummary.getKey().toString(),

                            summary.getSveInFile().getHashedMessage(), summary.getSveInFile().getReferenceAllele(),
                            summary.getSveInFile().getAlternateAllele(), summary.getSveInFile().getClusteredVariantAccession(),
                            summary.getRsInFileValid(), summary.getSveInFile().getRemappedFrom(),

                            summary.getSveInDB().getHashedMessage(), summary.getSveInDB().getReferenceAllele(),
                            summary.getSveInDB().getAlternateAllele(), summary.getSveInDB().getClusteredVariantAccession(),
                            summary.getRsInDBValid(), summary.getSveInDB().getRemappedFrom(),

                            summary.getCategory()}
                    br.write(row.join(","))
                    br.write("\n")
                }
            }
        }
    }

    void remediateOneValidRSWithRSInNonRemappedSVE() {
        logger.info("Processing SVEs where only one of them has a valid RS ")

        List<SubmittedVariantOperationEntity> svoeInsertList = new ArrayList<>()
        List<SubmittedVariantEntity> sveInsertList = new ArrayList<>()
        List<SubmittedVariantEntity> sveDeleteList = new ArrayList<>()

        List<SubmittedVariantOperationEntity> dbsnpSvoeInsertList = new ArrayList<>()
        List<SubmittedVariantEntity> dbsnpSveInsertList = new ArrayList<>()

        for (Map.Entry<String, List<CollisionSummary>> asmSummary : collisionSummaryMap.entrySet()) {
            String assembly = asmSummary.getKey()
            logger.info("processing for assembly : " + assembly)
            List<CollisionSummary> collisionSummaryList = asmSummary.getValue()
            // filter and take sve where only one valid rs is there
            List<CollisionSummary> impactedList = collisionSummaryList.stream()
                    .filter(colsum -> colsum.getCategory() == CATEGORY.BOTH_SVE_HAS_RS_ONE_VALID
                            || colsum.getCategory() == CATEGORY.ONE_SVE_HAS_RS_ONE_VALID)
                    .collect(Collectors.toList())
            logger.info("Impacted List Size : " + impactedList.size())

            // find our the residing collection
            Set<String> sveInFileHashes = impactedList.stream()
                    .map(colsum -> colsum.getSveInFile().getHashedMessage())
                    .collect(Collectors.toSet())
            Set<String> sveInDBHashes = impactedList.stream()
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

            for (CollisionSummary colSummary : impactedList) {
                //  There can be four different cases here
                //  1. if both are remapped or both are not remapped
                //      - we keep the one with the lower created date and sve accession
                //  2. if only one is remapped
                //      a. if the rs belongs to the non remapped keep that
                //      b. if the rs belongs to the remapped //TODO: confirm priority for the case

                SubmittedVariantEntity sveInFile = colSummary.getSveInFile()
                SubmittedVariantEntity sveInDB = colSummary.getSveInDB()
                ImmutablePair<Class, Class> sveInFileCollections = getVariantAndOperationCollection(sveInFile)
                ImmutablePair<Class, Class> sveInDBCollections = getVariantAndOperationCollection(sveInDB)

                // assume we need to keep the one in db
                Boolean sveToKeepIsInFile = false

                if ((sveInFile.getRemappedFrom() != null && sveInDB.getRemappedFrom() != null)
                        || (sveInFile.getRemappedFrom() == null && sveInDB.getRemappedFrom() == null)) {
                    if (sveInFile.getCreatedDate().isBefore(sveInDB.getCreatedDate())) {
                        sveToKeepIsInFile = true
                    } else if (sveInFile.getCreatedDate().isEqual(sveInDB.getCreatedDate())) {
                        if (sveInFile.getAccession() < sveInDB.getAccession()) {
                            sveToKeepIsInFile = true
                        }
                    }
                }

                if (sveToKeepIsInFile) {
                    // create a merge operation for sve in db
                    if (sveInDBCollections.getLeft() == sveClass) {
                        svoeInsertList.add(getSVOEForMergeOperation(sveInDB, sveInFile.getAccession(), false))
                    } else {
                        dbsnpSvoeInsertList.add(getSVOEForMergeOperation(sveInDB, sveInFile.getAccession(), false))
                    }
                    // delete the merged one from the db
                    sveDeleteList.add(sveInDB)
                    // remediate and add the sve in file to the db
                    if (sveInFileCollections.getLeft() == sveClass) {
                        sveInsertList.add(getRemediatedSVE(sveInFile))
                    } else {
                        dbsnpSveInsertList.add(getRemediatedSVE(sveInFile))
                    }

                } else {
                    // create a merge operation for sve in file
                    if (sveInFileCollections.getLeft() == sveClass) {
                        svoeInsertList.add(getSVOEForMergeOperation(sveInFile, sveInDB.getAccession(), true))
                    } else {
                        dbsnpSvoeInsertList.add(getSVOEForMergeOperation(sveInDB, sveInFile.getAccession(), false))
                    }

                    // delete the merged sve the db
                    sveDeleteList.add(sveInFile)
                }
            }
        }

        // add the merge operation for all SVEs
        logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Merged Operations List (SVOE): " + gson.toJson(svoeInsertList))
        logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Merged Operations List (DBSNPSVOE): " + gson.toJson(dbsnpSvoeInsertList))
        if (svoeInsertList != null && !svoeInsertList.isEmpty()) {
            // dbEnv.mongoTemplate.insert(svoeInsertList, svoeClass)
        }
        if (dbsnpSvoeInsertList != null && !dbsnpSvoeInsertList.isEmpty()) {
            // dbEnv.mongoTemplate.insert(dbsnpSvoeInsertList, dbsnpSvoeClass)
        }

        // delete the merged SVEs
        logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Delete Merged SVE : " + gson.toJson(sveDeleteList))
        // removeMergedSVE(sveDeleteList, sveClass)
        // removeMergedSVE(dbsnpSveDeleteList, dbsnpSveClass)

        // insert updated SVEs into DB
        logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Insert SVE : " + gson.toJson(sveInsertList))
        logger.info("remediateOneValidRSWithRSInNonRemappedSVE - Insert DBSNPSVE: " + gson.toJson(dbsnpSveInsertList))
        if (sveInsertList != null && !sveInsertList.isEmpty()) {
            // dbEnv.mongoTemplate.insert(sveInsertList, sveClass)
        }
        if (dbsnpSveInsertList != null && !dbsnpSveInsertList.isEmpty()) {
            // dbEnv.mongoTemplate.insert(dbsnpSveInsertList, dbsnpSveClass)
        }
    }

    ImmutablePair<Class, Class> getVariantAndOperationCollection(Set<String> sveCollection, Set<String> dbsnpSVECollection,
                                                                 SubmittedVariantEntity sve) {
        if (sveCollection.contains(sve.getHashedMessage())) {
            return new ImmutablePair(sveClass, svoeClass)
        } else if (dbsnpSVECollection.contains(sve.getHashedMessage())) {
            return new ImmutablePair(dbsnpSveClass, dbsnpSvoeClass)
        } else {
            throw new RuntimeException("Could not find the collection for sve " + gson.toJson(sve))
        }
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
        new SubmittedVariantEntity(submittedVariantEntity.getAccession(),
                getSVENewHash(submittedVariantEntity), sveModel, submittedVariantEntity.getVersion())
    }

    SubmittedVariantOperationEntity getSVOEForMergeOperation(SubmittedVariantEntity submittedVariantEntity, Long accessionMergeInto) {
        // remediate sve if it contains lowercase nucleotide before putting in svoe
        SubmittedVariantEntity sve = getRemediatedSVE(submittedVariantEntity)

        // create merge operation
        SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity()
        svoe.fill(EventType.MERGED, sve.getAccession(), accessionMergeInto,
                "After fixing lowercase nucleotide issue, variant merged due to duplicate hash",
                Collections.singletonList(new SubmittedVariantInactiveEntity(submittedVariantEntity)))
        svoe.setId("EVA3417_MERGED_" + sve.getAccession() + "_INTO_" + accessionMergeInto + "_HASH_" + sve.getHashedMessage())

        return svoe
    }

    void runSummary() {
        buildSummary()
        //printSummary()
        writeSummaryToFile()

        // remediate the case where only one sve has a valid rs
        remediateOneValidRSWithRSInNonRemappedSVE()
    }
}


enum CATEGORY {
    NO_SVE_HAS_RS,
    ONE_SVE_HAS_RS_NONE_VALID,
    ONE_SVE_HAS_RS_ONE_VALID,  // around 23262
    BOTH_SVE_HAS_RS_NONE_VALID,
    BOTH_SVE_HAS_RS_ONE_VALID,  // around 17
    BOTH_SVE_HAS_RS_BOTH_VALID  // around 221
}

class CollisionSummary {
    private SubmittedVariantEntity sveInFile
    private Long rsInFile
    private boolean rsInFileValid
    private List<ClusteredVariantEntity> fileRSList

    private SubmittedVariantEntity sveInDB
    private Long rsInDB
    private boolean rsInDBValid
    private List<ClusteredVariantEntity> dbRSList

    private CATEGORY category

    CollisionSummary(SubmittedVariantEntity sveInFile, SubmittedVariantEntity sveInDB, Map<Long, List<ClusteredVariantEntity>> cveInDB) {
        this.sveInFile = sveInFile
        this.rsInFile = sveInFile.getClusteredVariantAccession()
        this.fileRSList = cveInDB.get(rsInFile)
        if (fileRSList == null || fileRSList.isEmpty()) {
            this.rsInFileValid = false
        } else {
            this.rsInFileValid = true
        }

        this.sveInDB = sveInDB
        this.rsInDB = sveInDB.getClusteredVariantAccession()
        this.dbRSList = cveInDB.get(rsInDB)
        if (dbRSList == null || dbRSList.isEmpty()) {
            this.rsInDBValid = false
        } else {
            this.rsInDBValid = true
        }

        this.category = identifyCategory()
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

    Long getRsInFile() {
        return rsInFile
    }

    boolean getRsInFileValid() {
        return rsInFileValid
    }

    List<ClusteredVariantEntity> getFileRSList() {
        return fileRSList
    }

    SubmittedVariantEntity getSveInDB() {
        return sveInDB
    }

    Long getRsInDB() {
        return rsInDB
    }

    boolean getRsInDBValid() {
        return rsInDBValid
    }

    List<ClusteredVariantEntity> getDbRSList() {
        return dbRSList
    }

    CATEGORY getCategory() {
        return category
    }
}


