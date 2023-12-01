package eva3417

import com.google.gson.Gson
import com.google.gson.reflect.TypeToken
import groovy.cli.picocli.CliBuilder
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Criteria
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.summary.SubmittedVariantSummaryFunction
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment

import java.util.function.Function
import java.util.stream.Collectors

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*


def cli = new CliBuilder()
cli.envPropertiesFile(args: 1, "properties file for connecting to db", required: true)
cli.logFile(args: 2, "Path to the log file for running QC", required: true)

def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}


// call to run QC
new RemediationQC(options.envPropertiesFile, options.logFile).runQC()

class RemediationQC {
    private String logFile
    private EVADatabaseEnvironment dbEnv
    private Function<ISubmittedVariant, String> hashingFunction
    static def logger = LoggerFactory.getLogger(Application.class)
    private Gson gson

    RemediationQC(envPropertiesFile, logFile) {
        this.dbEnv = createFromSpringContext(envPropertiesFile, Application.class)
        this.logFile = logFile
        this.hashingFunction = new SubmittedVariantSummaryFunction().andThen(new SHA1HashingFunction())
        this.gson = new Gson()
    }

    String getSVENewHash(SubmittedVariantEntity sve) {
        SubmittedVariant temp = sve.getModel()
        temp.setReferenceAllele(sve.getReferenceAllele().toUpperCase())
        temp.setAlternateAllele(sve.getAlternateAllele().toUpperCase())
        return hashingFunction.apply(temp)
    }


    ArrayList<SubmittedVariantEntity> getListOfSVEFromLine(String linePrefix, String line) {
        String sveListString = line.substring(line.indexOf(linePrefix) + linePrefix.length())
        ArrayList<SubmittedVariantEntity> sveList = (ArrayList<SubmittedVariantEntity>) gson.fromJson(sveListString,
                new TypeToken<ArrayList<SubmittedVariantEntity>>() {}.getType());

        return sveList
    }

    ArrayList<SubmittedVariantOperationEntity> getListOfSVOEFromLine(String linePrefix, String line){
        String svoeListString = line.substring(line.indexOf(linePrefix) + linePrefix.length())
        ArrayList<SubmittedVariantOperationEntity> svoeList = (ArrayList< SubmittedVariantOperationEntity>) gson.fromJson(svoeListString,
                new TypeToken<ArrayList<SubmittedVariantOperationEntity>>() {}.getType())
        return svoeList
    }

    void checkAllSVERemediatedAndInserted(List<SubmittedVariantEntity> sveList, String logPrefix){
        Set<String> newHashes = sveList.stream().map(sve -> getSVENewHash(sve))
                .collect(Collectors.toSet())
        List<SubmittedVariantEntity> newHashSVEInDB = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(newHashes)), ssClass))
                .flatten()
        // check if all the SVEs inserted are remediated
        List<SubmittedVariantEntity> notRemediatedList = newHashSVEInDB.stream()
                .filter(sve -> sve.getReferenceAllele() != sve.getReferenceAllele().toUpperCase()
                        || sve.getAlternateAllele() != sve.getAlternateAllele().toUpperCase())
                .collect(Collectors.toList())
        if (notRemediatedList.size() > 0) {
            logger.error(logPrefix + " (SVE not remediated): " + gson.toJson(notRemediatedList))
        }
        // check if all the SVE were inserted
        if (newHashSVEInDB.size() < sveList.size()) {
            Set<String> hashesInDB = newHashSVEInDB.stream().map(sve -> sve.getHashedMessage())
                    .collect(Collectors.toSet())
            List<SubmittedVariantEntity> missingEntries = new ArrayList<>()
            for (SubmittedVariantEntity sve : sveList) {
                if (!hashesInDB.contains(getSVENewHash(sve))) {
                    missingEntries.add(sve)
                }
            }
            logger.error(logPrefix + " (Entries Missing for SVE): " + gson.toJson(missingEntries))
        }
    }

    void checkAllSVEWereDeletedFromDB(List<SubmittedVariantEntity> sveList, String logPrefix){
        Set<String> oldHashes = sveList.stream().map(sve -> sve.getHashedMessage())
                .collect(Collectors.toSet())
        List<SubmittedVariantEntity> oldHashSVEList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(oldHashes)), ssClass))
                .flatten()
        if(oldHashSVEList.size() > 0){
            logger.error(logPrefix + gson.toJson(oldHashSVEList))
        }
    }

    void qcNoHashCollisionEntries(List<SubmittedVariantEntity> sveList) {
        // check if all SVEs were remediated and inserted
        checkAllSVERemediatedAndInserted(sveList, "No Hash Collision")

        // check if all the existing SVE were deleted after remediation
        checkAllSVEWereDeletedFromDB(sveList, "No Hash Collision (SVE not deleted after remediation) : ")

        // check if all the update operations were created
        Set<String> updateOpIds = sveList.stream()
                .map(sve -> "EVA3417_UPDATED_" + sve.getReferenceSequenceAccession() + "_" + sve.getAccession())
                .collect(Collectors.toSet())
        List<SubmittedVariantOperationEntity> updateOpInDB = [svoeClass, dbsnpSvoeClass]
                .collectMany(opClass -> dbEnv.mongoTemplate.find(query(where("_id").in(updateOpIds)), opClass))
                .flatten()
        if (updateOpInDB.size() != sveList.size()) {
            Set<String> updateOPIdsInDB = updateOpInDB.stream().map(updateOp -> updateOp.getId())
                    .collect(Collectors.toSet())
            List<SubmittedVariantOperationEntity> missingEntries = new ArrayList<>()
            for (SubmittedVariantEntity sve : sveList) {
                if (!updateOPIdsInDB.contains("EVA3417_UPDATED_" + sve.getReferenceSequenceAccession() + "_" + sve.getAccession())) {
                    missingEntries.add(sve)
                }
            }
            logger.error("No Hash Collision (Entries Missing in SVOE for update Operation): " + gson.toJson(missingEntries))
        }
    }

    void qcHashCollisionDiscardEntries(List<SubmittedVariantEntity> sveList) {
        // check none of the discarded SVEs are present in DB
        checkAllSVEWereDeletedFromDB(sveList, "Hash Collision (discarded SVE not deleted): ")

        // check if all the discard operations were created
        Set<String> discardOpIds = sveList.stream()
                .map(sve -> "EVA3417_DISCARD_SS_" + sve.getAccession() + "_HASH_" + sve.getHashedMessage())
                .collect(Collectors.toSet())
        List<SubmittedVariantOperationEntity> discardOpInDB = [svoeClass, dbsnpSvoeClass]
                .collectMany(opClass -> dbEnv.mongoTemplate.find(query(where("_id").in(discardOpIds)), opClass))
                .flatten()
        if (discardOpInDB.size() != sveList.size()) {
            Set<String> discardOPIdsInDB = discardOpInDB.stream().map(discardOp -> discardOp.getId())
                    .collect(Collectors.toSet())
            List<SubmittedVariantOperationEntity> missingEntries = new ArrayList<>()
            for (SubmittedVariantEntity sve : sveList) {
                if (!discardOPIdsInDB.contains("EVA3417_DISCARD_SS_" + sve.getAccession() + "_HASH_" + sve.getHashedMessage())) {
                    missingEntries.add(sve)
                }
            }
            logger.error("Hash Collision (Entries Missing in SVOE for discard operation): " + gson.toJson(missingEntries))
        }
    }

    void qcHashCollisionInsertEntries(List<SubmittedVariantEntity> sveList){
        checkAllSVERemediatedAndInserted(sveList, "Hash Collision")
    }

    void qcHashCollisionRemoveEntries(List<SubmittedVariantEntity> sveList){
        // check if all the SVEs were removed
        checkAllSVEWereDeletedFromDB(sveList, "Hash Collision (SVE not removed): ")
    }

    void qcHashCollisionMergeEntries(List<SubmittedVariantEntity> sveList){
        // check if all the merged SVE were deleted
        Set<String> oldHashes = sveList.stream().map(sve -> sve.getHashedMessage())
                .collect(Collectors.toSet())
        Set<Long> oldAccessions = sveList.stream().map(sve->sve.getAccession())
                .collect(Collectors.toSet())
        List<SubmittedVariantEntity> oldHashSVEList = [sveClass, dbsnpSveClass]
                .collectMany(ssClass -> dbEnv.mongoTemplate.find(query(where("_id").in(oldHashes).and("accession").in(oldAccessions)), ssClass))
                .flatten()
        if (oldHashSVEList.size() > 0) {
            logger.error("Hash Collision (SVE not deleted after merge): " + gson.toJson(oldHashSVEList))
        }

        // check if all the merge operations were created
        Set<String> mergedHashes =  sveList.stream()
                .map(sve->getSVENewHash(sve))
                .collect(Collectors.toSet())
        Set<String> mergedAssemblies = sveList.stream()
                .map(sve->sve.getReferenceSequenceAccession())
                .collect(Collectors.toSet())
        Set<Long> mergedAccessions = sveList.stream()
                .map(sve->sve.getAccession())
                .collect(Collectors.toSet())
        Criteria filterCriteria = new Criteria("inactiveObjects.seq").in(mergedAssemblies)
                .and("accession").in(mergedAccessions)
                .and("eventType").is("MERGED")
                .and("inactiveObjects.hashedMessage").in(mergedHashes)
                .and("_id").regex("^EVA3417_MERGED_")

        List<SubmittedVariantOperationEntity> mergeOpInDB = [svoeClass, dbsnpSvoeClass]
                        .collectMany(opClass -> dbEnv.mongoTemplate.find(query(filterCriteria), opClass))
                        .flatten()
        if (mergeOpInDB.size() != sveList.size()) {
            //TODO: filter out and log only the missing entries
            logger.error("Hash Collision (Entries Missing in SVOE for Merge Operation): " + gson.toJson(sveList))
        }
    }

    void qcExistingSVOEOperationAreRemediated(List<SubmittedVariantOperationEntity> svoeList){
        List<String> svoeIds = svoeList.stream().map(svoe->svoe.getId())
                .collect(Collectors.toList())
        List<SubmittedVariantOperationEntity> svoeInDB = [svoeClass, dbsnpSvoeClass]
                .collectMany(opClass -> dbEnv.mongoTemplate.find(query(where("_id").in(svoeIds)), opClass))
                .flatten()
        List<SubmittedVariantOperationEntity> svoeNotUpdated = new ArrayList<>();
        for(SubmittedVariantOperationEntity svoe: svoeInDB){
            SubmittedVariantInactiveEntity inactiveSVEDocument = svoe.getInactiveObjects()[0]
            if(inactiveSVEDocument.getReferenceAllele()!=inactiveSVEDocument.getReferenceAllele().toUpperCase()
                    ||inactiveSVEDocument.getAlternateAllele()!=inactiveSVEDocument.getAlternateAllele()){
                svoeNotUpdated.add(svoe)
            }
        }

        if(svoeNotUpdated.size()>0){
            logger.error("Existing SVOE not updated: " + gson.toJson(svoeNotUpdated))
        }
    }

    void runQC() {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(this.logFile)))) {
            String line
            while ((line = br.readLine()) != null) {
                if (line.contains("Impacted sve List (No Hash Collision): ")) {
                    List<SubmittedVariantEntity> sveList = getListOfSVEFromLine("Impacted sve List (No Hash Collision): ", line)
                    qcNoHashCollisionEntries(sveList)
                } else if (line.contains("Hash Collision Discard SVE List: ")) {
                    List<SubmittedVariantEntity> sveList = getListOfSVEFromLine("Hash Collision Discard SVE List: ", line)
                    qcHashCollisionDiscardEntries(sveList)
                } else if (line.contains("Hash Collision Merge SVE List: ")) {
                    List<SubmittedVariantEntity> sveList = getListOfSVEFromLine("Hash Collision Merge SVE List: ", line)
                    qcHashCollisionMergeEntries(sveList)
                } else if (line.contains("Hash Collision Insert SVE List: ")) {
                    List<SubmittedVariantEntity> sveList = getListOfSVEFromLine("Hash Collision Insert SVE List: ", line)
                    qcHashCollisionInsertEntries(sveList)
                } else if (line.contains("Hash Collision Remove SVE List: ")) {
                    List<SubmittedVariantEntity> sveList = getListOfSVEFromLine("Hash Collision Remove SVE List: ", line)
                    qcHashCollisionRemoveEntries(sveList)
                } else if (line.contains("Existing SVOE operations to update: ")) {
                    List<SubmittedVariantOperationEntity> svoeList = getListOfSVOEFromLine("Existing SVOE operations to update: ", line)
                    qcExistingSVOEOperationAreRemediated(svoeList)
                }
            }
        } catch (IOException e) {
            // Handle the exception appropriately
            e.printStackTrace()
        }
    }

}
