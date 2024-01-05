package eva3417

import com.google.gson.Gson
import com.google.gson.reflect.TypeToken
import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory
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

new SummariseHashCollision(options.envPropertiesFile, options.baseDirPath, options.summaryFilePath).runSummaryAndRemediate()

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
                for (int i = 0; i < sveList.size(); i++) {
                    SubmittedVariantEntity sveInFile = sveList.get(i)
                    SubmittedVariantEntity sveInDB = sveInDBMap.get(getSVENewHash(sveInFile))
                    collisionSummaryMap.get(asm).add(new CollisionSummary(sveInFile, sveInDB, cveInDB))
                }
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

    SubmittedVariantEntity getSveForHash(String hash) {
        List<SubmittedVariantEntity> sveList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").is(hash)), ssClass))
                .flatten()
        if (sveList.size() > 1) {
            throw new RuntimeException("Multiple sve found for hash " + hash)
        }
        return sveList.get(0)
    }

    Class getCollectionOfSVE(String hash) {
        List<SubmittedVariantEntity> sveList = dbEnv.mongoTemplate.find(query(where("_id").is(hash)), sveClass)
        if (sveList != null && !sveList.isEmpty()) {
            return sveClass
        } else {
            List<SubmittedVariantEntity> dbsnpSveList = dbEnv.mongoTemplate.find(query(where("_id").is(hash)), dbsnpSveClass)
            if (dbsnpSveList != null && !dbsnpSveList.isEmpty()) {
                return dbsnpSveClass
            } else {
                throw new RuntimeException("SVE with hash " + hash + " could not be found")
            }
        }
    }

    List<ClusteredVariantEntity> getCVEInDB(Long... accessions) {
        Set<Long> accList = new ArrayList<>()
        for (Long acc : accessions) {
            accList.add(acc)
        }
        logger.info("getting CVE for accessions " + gson.toJson(accList))
        List<ClusteredVariantEntity> cveInDBList = [cveClass, dbsnpCveClass]
                .collectMany(cveClass -> dbEnv.mongoTemplate.find(query(where("accession").in(accList)), cveClass))
                .flatten()

        return cveInDBList
    }

    boolean checkIfValidRS(Long rs, String asm, List<ClusteredVariantEntity> cveList) {
        if (rs == null) {
            return false
        } else {
            Optional<ClusteredVariantEntity> opRS = cveList.stream()
                    .filter(cve -> cve.getAccession() == rs && cve.getAssemblyAccession() == asm)
                    .findAny()
            if (opRS.isPresent()) {
                return true
            } else {
                return false
            }
        }
    }

    static CATEGORY getCollisionCategory(Long rsInFile, Long rsInDB, boolean rsInFileValid, boolean rsInDBValid) {
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

    List<SubmittedVariantEntity> findSVEAccessionInOrigAssembly(SubmittedVariantEntity sve) {
        List<SubmittedVariantEntity> dbSVEList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("seq").is(sve.getRemappedFrom()).and("accession").is(sve.getAccession())), ssClass))
                .flatten()
        return dbSVEList
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


    void remediateCategoryAtMostOneValidRS() {
        String[] affectedAsmList = [
                "GCA_000002285.2", "GCA_000002315.3", "GCA_000003055.5", "GCA_011100555.1", "GCA_000002035.3",
                "GCA_000001635.6", "GCA_000349105.1", "GCA_016772045.1", "GCA_000003025.6", "GCA_000349185.1",
                "GCA_000247795.2", "GCA_000473445.2",
                "GCA_001433935.1"
        ]
        for (String assembly : affectedAsmList) {
            logger.info("processing for assembly: " + assembly)
            File asmHashCollisionFile = new File(baseDirPath + "/" + assembly + "/hash_collision.txt")
            // read all the hash collision SVE from the file
            List<SubmittedVariantEntity> sveInCollisionFile = Files.lines(asmHashCollisionFile.toPath())
                    .flatMap(line -> getListOfSVEFromLine(line).stream())
                    .collect(Collectors.toList())
            logger.info("Number of sve in collision file: " + sveInCollisionFile.size())
            // process all SVEs one by one
            int varNum = 1
            for (SubmittedVariantEntity sveInFile : sveInCollisionFile) {
                logger.info("processing variant number: " + varNum)
                varNum++
                logger.info("processing sve in file: " + gson.toJson(sveInFile))
                SubmittedVariantEntity sveInDB = getSveForHash(getSVENewHash(sveInFile))
                logger.info("processing sve in DB: " + gson.toJson(sveInDB))

                Class sveInFileCollection = getCollectionOfSVE(sveInFile.getHashedMessage())
                Class sveInDBCollection = getCollectionOfSVE(sveInDB.getHashedMessage())
                logger.info("sveInFileCollection: " + sveInFileCollection)
                logger.info("sveInDBCollection: " + sveInDBCollection)

                Long rsInFile = sveInFile.getClusteredVariantAccession()
                Long rsInDB = sveInDB.getClusteredVariantAccession()
                logger.info("rsInFile: " + rsInFile)
                logger.info("rsInDB: " + rsInDB)

                List<ClusteredVariantEntity> cveInDBList = getCVEInDB(rsInFile, rsInDB)
                logger.info("CVE list: " + gson.toJson(cveInDBList))

                boolean rsInFileValid = checkIfValidRS(rsInFile, sveInFile.getReferenceSequenceAccession(), cveInDBList)
                boolean rsInDBValid = checkIfValidRS(rsInDB, sveInDB.getReferenceSequenceAccession(), cveInDBList)
                logger.info("rsInFile Valid: " + rsInFileValid)
                logger.info("rsInDB Valid: " + rsInDBValid)

                CATEGORY colCategory = getCollisionCategory(rsInFile, rsInDB, rsInFileValid, rsInDBValid)
                logger.info("Collision Category: " + colCategory)

                // check if only one valid rs is present among both collision SVEs
                if (colCategory == CATEGORY.BOTH_SVE_HAS_RS_BOTH_VALID) {
                    logger.info("Collision Category is having both RS Valid. skipping")
                    continue
                }

                // assume we need to keep the sve in db and merge the sve in file (lowercase)
                SubmittedVariantEntity sveToKeep = sveInDB
                SubmittedVariantEntity sveToMerge = sveInFile
                boolean sveToKeepRSValid = rsInDBValid
                boolean sveToMergeRSValid = rsInFileValid
                Boolean sveToKeepIsInFile = false

                if (sveInFile.getRemappedFrom() == null && sveInDB.getRemappedFrom() == null) {
                    // case 1: both are not remapped (take the one with lowest created date/accession)
                    if (sveInFile.getCreatedDate().isBefore(sveInDB.getCreatedDate())) {
                        sveToKeep = sveInFile
                        sveToMerge = sveInDB
                        sveToKeepRSValid = rsInFileValid
                        sveToMergeRSValid = rsInDBValid
                        sveToKeepIsInFile = true
                    } else if (sveInFile.getCreatedDate().isEqual(sveInDB.getCreatedDate())) {
                        if (sveInFile.getAccession() < sveInDB.getAccession()) {
                            sveToKeep = sveInFile
                            sveToMerge = sveInDB
                            sveToKeepRSValid = rsInFileValid
                            sveToMergeRSValid = rsInDBValid
                            sveToKeepIsInFile = true
                        }
                    }
                } else if (sveInFile.getRemappedFrom() != null && sveInDB.getRemappedFrom() != null) {
                    // case 2: both are remapped check which sve accession is present in original assembly that gets the priority
                    List<SubmittedVariantEntity> sveInFileOrgAsm = findSVEAccessionInOrigAssembly(sveInFile)
                    List<SubmittedVariantEntity> sveInDBOrgAsm = findSVEAccessionInOrigAssembly(sveInDB)
                    if (((sveInFileOrgAsm != null && !sveInFileOrgAsm.isEmpty()) && (sveInDBOrgAsm != null && !sveInDBOrgAsm.isEmpty()))
                            || ((sveInFileOrgAsm == null || sveInFileOrgAsm.isEmpty()) && (sveInDBOrgAsm == null || sveInDBOrgAsm.isEmpty()))) {
                        // both sve accession are present in original assembly or both sve accession are not present in the original assembly
                        // similar to case 1 - take the one with lower created date and sve accession
                        if (sveInFile.getCreatedDate().isBefore(sveInDB.getCreatedDate())) {
                            sveToKeep = sveInFile
                            sveToMerge = sveInDB
                            sveToKeepRSValid = rsInFileValid
                            sveToMergeRSValid = rsInDBValid
                            sveToKeepIsInFile = true
                        } else if (sveInFile.getCreatedDate().isEqual(sveInDB.getCreatedDate())) {
                            if (sveInFile.getAccession() < sveInDB.getAccession()) {
                                sveToKeep = sveInFile
                                sveToMerge = sveInDB
                                sveToKeepRSValid = rsInFileValid
                                sveToMergeRSValid = rsInDBValid
                                sveToKeepIsInFile = true
                            }
                        }
                    } else {
                        // only one sve accession is present in original assembly, take the one that is present
                        if ((sveInFileOrgAsm != null && !sveInFileOrgAsm.isEmpty()) && (sveInDBOrgAsm == null || sveInDBOrgAsm.isEmpty())) {
                            sveToKeep = sveInFile
                            sveToMerge = sveInDB
                            sveToKeepRSValid = rsInFileValid
                            sveToMergeRSValid = rsInDBValid
                            sveToKeepIsInFile = true
                        }
                    }
                } else if (sveInFile.getRemappedFrom() == null || sveInDB.getRemappedFrom() == null) {
                    // case 3: only one is remapped , non remapped one gets the priority
                    if (sveInFile.getRemappedFrom() == null && sveInDB.getRemappedFrom() != null) {
                        sveToKeep = sveInFile
                        sveToMerge = sveInDB
                        sveToKeepRSValid = rsInFileValid
                        sveToMergeRSValid = rsInDBValid
                        sveToKeepIsInFile = true
                    }
                }

                logger.info("sve to keep: " + gson.toJson(sveToKeep))
                logger.info("sve to merge: " + gson.toJson(sveToMerge))

                Class sveToKeepCollection = getCollectionOfSVE(sveToKeep.getHashedMessage())
                Class sveToMergeCollection = getCollectionOfSVE(sveToMerge.getHashedMessage())
                logger.info("sve to keep collection: " + sveToKeepCollection)
                logger.info("sve to merge collection: " + sveToMergeCollection)

                // create a merge operation for SVE to merge
                SubmittedVariantOperationEntity sveMergeOperation = getSVOEForMergeOperation(sveToMerge, sveToKeep.getAccession())
                logger.info("SVOE Merge Collection Entry : " + gson.toJson(sveMergeOperation))
                dbEnv.mongoTemplate.insert(Arrays.asList(sveMergeOperation), sveToMergeCollection == sveClass ? svoeClass : dbsnpSvoeClass)

                // delete the merged SVE
                logger.info("Delete Merged SVE : " + gson.toJson(sveToMerge))
                Query mergedDeleteQuery = query(where("_id").is(sveToMerge.getHashedMessage()).and("accession").in(sveToMerge.getAccession()))
                dbEnv.mongoTemplate.remove(mergedDeleteQuery, sveToMergeCollection)

                if (sveToKeepIsInFile) {
                    // delete the lowercase sve (will insert the remediated one)
                    logger.info("Delete LowerCase SVE: " + gson.toJson(sveToKeep))
                    Query lowerCaseDeleteQuery = query(where("_id").is(sveToKeep.getHashedMessage()).and("accession").in(sveToKeep.getAccession()))
                    dbEnv.mongoTemplate.remove(lowerCaseDeleteQuery, sveToKeepCollection)

                    // if SVE to keep does not have a valid RS but SVE to merge has, copy the valid RS to SVE to Keep
                    if (sveToKeepRSValid == false && sveToMergeRSValid == true) {
                        sveToKeep.setClusteredVariantAccession(sveToMerge.getClusteredVariantAccession())
                    }

                    // remediate and add the sve to the db
                    SubmittedVariantEntity remediatedSVE = getRemediatedSVE(sveToKeep)
                    logger.info("Remediated SVE: " + gson.toJson(remediatedSVE))
                    dbEnv.mongoTemplate.insert(Arrays.asList(remediatedSVE), sveToKeepCollection)
                } else {
                    // if sve to keep is present in db but does not have a valid rs, then update the rs
                    if (sveToKeepRSValid == false && sveToMergeRSValid == true) {
                        logger.info("SVE for which RS needs to be updated: " + gson.toJson(sveToKeep))
                        sveToKeep.setClusteredVariantAccession(sveToMerge.getClusteredVariantAccession())
                        logger.info("SVE after RS update: " + gson.toJson(sveToKeep))

                        Query updateFindQuery = query(where("_id").is(sveToKeep.getHashedMessage()))
                        Update updateQuery = new Update().set("rs", sveToKeep.getClusteredVariantAccession())
                        dbEnv.mongoTemplate.updateFirst(updateFindQuery, updateQuery, sveToKeepCollection)
                    }

                }
            }
        }
    }

    void runSummaryAndRemediate() {
        buildSummary()
        //printSummary()
        writeSummaryToFile()

        // remediate the case where involved sve has at most one valid rs
        remediateCategoryAtMostOneValidRS()
    }
}


enum CATEGORY {
    NO_SVE_HAS_RS,
    ONE_SVE_HAS_RS_NONE_VALID,
    ONE_SVE_HAS_RS_ONE_VALID,
    BOTH_SVE_HAS_RS_NONE_VALID,
    BOTH_SVE_HAS_RS_ONE_VALID,
    BOTH_SVE_HAS_RS_BOTH_VALID
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

        this.category = SummariseHashCollision.getCollisionCategory(rsInFile, rsInDB, rsInFileValid, rsInDBValid)
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


