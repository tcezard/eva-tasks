package eva3372

import org.slf4j.LoggerFactory
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.getDbsnpSveClass
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.getSveClass

class CollectIdsAndHashesToFile {
    static def logger = LoggerFactory.getLogger(CollectIdsAndHashesToFile.class)

    String assembly
    EVADatabaseEnvironment prodEnv
    String idAndHashFilePath

    CollectIdsAndHashesToFile() {}

    CollectIdsAndHashesToFile(String assembly, EVADatabaseEnvironment prodEnv, String idAndHashFilePath) {
        this.assembly = assembly
        this.prodEnv = prodEnv
        this.idAndHashFilePath = idAndHashFilePath
    }

    def exportRsIdsAndHashes = {
        def evaAndDbsnpSveCursorsProd = [sveClass, dbsnpSveClass].collect { collectionClass ->
            new RetryableBatchingCursor<>(
                    where("seq").is(assembly).and("rs").exists(true),
                    prodEnv.mongoTemplate, collectionClass)
        }
        def numRecordsProcessed = 0
        evaAndDbsnpSveCursorsProd.each { cursor ->
            cursor.each { List<SubmittedVariantEntity> sves ->
                def remappedSves = sves.findAll { !it.remappedFrom.equals("") }
                if (!remappedSves.isEmpty()) {
                    def rsIdAndCveHash = remappedSves.collect{ new Tuple(it.clusteredVariantAccession, EVAObjectModelUtils.getClusteredVariantHash(it)) }
                    writeToFile(rsIdAndCveHash)
                    numRecordsProcessed += remappedSves.size()
                    logger.info("${numRecordsProcessed} SS processed so far...")
                }
            }
        }
        logger.info("Scan through SVE and dbsnpSVE collections complete")
    }

    def writeToFile = { List<Tuple> rsIdAndCveHash ->
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(idAndHashFilePath, true))) {
            rsIdAndCveHash.each { writer.writeLine("${it[0]}\t${it[1]}") }
        } catch (IOException e) {
            e.printStackTrace()
        }
    }

}
