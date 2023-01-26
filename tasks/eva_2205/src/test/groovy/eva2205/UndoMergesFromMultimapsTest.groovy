package eva2205

import org.junit.After
import org.junit.Before
import org.junit.Test
import org.junit.runner.RunWith
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.context.ApplicationContext
import org.springframework.test.context.TestPropertySource
import org.springframework.test.context.junit4.SpringRunner
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType

import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.CommonTestUtils
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.EVAObjectModelUtils

import static org.junit.Assert.assertEquals
import static org.junit.Assert.assertTrue
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where

@RunWith(SpringRunner.class)
@TestPropertySource("classpath:application-test.properties")
class UndoMergesFromMultimapsTest {
    @Autowired
    private ApplicationContext applicationContext

    private boolean dbEnvSetupDone = false

    private static final String ASSEMBLY = "GCA_000000001.1"

    private static final int TAXONOMY = 60711

    private static EVADatabaseEnvironment dbEnv

    private static List<SubmittedVariantOperationEntity> svoeMerges

    private static void cleanup() {
        dbEnv.mongoClient.dropDatabase(dbEnv.mongoTemplate.db.name)
        svoeMerges = new ArrayList<SubmittedVariantOperationEntity>()
    }

    @Before
    void setUp() {
        if (!this.dbEnvSetupDone) {
            // We need to directly use the application-test.properties file here
            // because there are symbolic variables like |eva.mongo.host.test|
            // that needs to be resolved by Spring and Maven dynamically
            def dynamicPropertyFilePath = CommonTestUtils.getTempPropertyFilePath(applicationContext)
            dbEnv = createFromSpringContext(dynamicPropertyFilePath, GenericApplication.class)
            this.dbEnvSetupDone = true
        }
        cleanup()
    }

    @After
    void tearDown() {
        cleanup()
    }

    private static def createSS(Long ssID, Long rsID, Map rsLocus, Integer mapWeight = null) {
        def svObj = new SubmittedVariant(ASSEMBLY, TAXONOMY, "prj1", "chr1",
                rsLocus.start, rsLocus.ref, rsLocus.alt, rsID)
        if (Objects.nonNull(mapWeight)) svObj.setMapWeight(mapWeight)
        return EVAObjectModelUtils.toSubmittedVariantEntity(ssID, svObj)
    }

    private static def mergeSvoeOp(SubmittedVariantEntity impactedSS, Long destRSID) {
        def op = new SubmittedVariantOperationEntity()
        op.fill(EventType.UPDATED, impactedSS.accession, "Original rs${impactedSS.clusteredVariantAccession} " +
                "was merged into rs${destRSID}.", Arrays.asList(new SubmittedVariantInactiveEntity(impactedSS)))
        return op
    }

    private void setupScenario1() {
        def rsLocus1  = ["start": 100L, "ref": "A", "alt": "C"]
        def rsLocus2  = ["start": 101L, "ref": "T", "alt": "G"]
        def ss1 = createSS(1, 1, rsLocus1)
        def ss2 = createSS(2, 1, rsLocus2)
        def ss2WithRS2 = createSS(2, 2, rsLocus2, 3)
        def ss2WithRS3 = createSS(2, 3, rsLocus2, 3)
        def (rs1WithLocus1, rs1WithLocus2) = [ss1, ss2].collect {
            EVAObjectModelUtils.toClusteredVariantEntity(it) }

        def svoeOp1 = mergeSvoeOp(ss2WithRS3, 2)
        def svoeOp2 = mergeSvoeOp(ss2WithRS2, 1)

        dbEnv.bulkInsertIgnoreDuplicates([ss1, ss2], dbsnpSveClass)
        dbEnv.bulkInsertIgnoreDuplicates([rs1WithLocus1, rs1WithLocus2], dbsnpCveClass)
        dbEnv.bulkInsertIgnoreDuplicates([svoeOp1, svoeOp2], dbsnpSvoeClass)
        // To simplify, we don't create any dbsnpCvoe records since it is not used by the undo script
    }

    private void setupScenario2() {
        def rsLocus1  = ["start": 100L, "ref": "A", "alt": "C"]
        def rsLocus2  = ["start": 101L, "ref": "T", "alt": "G"]
        def ss1 = createSS(1, 1, rsLocus1)
        def ss2 = createSS(2, 1, rsLocus2)
        def ss2WithRS2 = createSS(2, 2, rsLocus2, 3)
        def ss2WithRS3 = createSS(2, 3, rsLocus2)
        def (rs1WithLocus1, rs1WithLocus2) = [ss1, ss2].collect {
            EVAObjectModelUtils.toClusteredVariantEntity(it) }

        def svoeOp1 = mergeSvoeOp(ss2WithRS3, 2)
        def svoeOp2 = mergeSvoeOp(ss2WithRS2, 1)

        dbEnv.bulkInsertIgnoreDuplicates([ss1, ss2], dbsnpSveClass)
        dbEnv.bulkInsertIgnoreDuplicates([rs1WithLocus1, rs1WithLocus2], dbsnpCveClass)
        dbEnv.bulkInsertIgnoreDuplicates([svoeOp1, svoeOp2], dbsnpSvoeClass)
        // To simplify, we don't create any dbsnpCvoe records since it is not used by the undo script
    }

    @Test
    // See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=1591315552
    void testMergeFromMultimapScenario1() {
        setupScenario1()
        new UndoMergesFromMultiMaps(dbEnv, ASSEMBLY).undoMergesFromMultimapRS()

        // Ensure ss2 is updated with rs2 and mapWeight 3
        def updatedRecordForSS2 = dbEnv.mongoTemplate.find(query(where("accession").is(2L)), dbsnpSveClass)[0]
        assertEquals(3, updatedRecordForSS2.mapWeight)
        assertEquals(2L, updatedRecordForSS2.clusteredVariantAccession)

        // Ensure the action of SS2 being assigned rs2 above is recorded
        def opRecordForSS2MapWt = dbEnv.mongoTemplate.find(query(where("inactiveObjects.hashedMessage")
                .is(updatedRecordForSS2.hashedMessage).and("_id").regex("RESTORE_MAPWT_.*")), dbsnpSvoeClass)[0]
        assertEquals("RESTORE_MAPWT_" + updatedRecordForSS2.hashedMessage, opRecordForSS2MapWt.id)
        assertEquals(1L, opRecordForSS2MapWt.inactiveObjects[0].clusteredVariantAccession)
        assertTrue(Objects.isNull(opRecordForSS2MapWt.inactiveObjects[0].mapWeight))
    }

    @Test
    // See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=1898709629
    void testMergeFromMultimapScenario2() {
        setupScenario2()
        new UndoMergesFromMultiMaps(dbEnv, ASSEMBLY).undoMergesFromMultimapRS()

        // Ensure ss2 still remains
        def recordForSS2 = dbEnv.mongoTemplate.find(query(where("accession").is(2L)), dbsnpSveClass)[0]
        assertTrue(Objects.isNull(recordForSS2.mapWeight))
        assertEquals(1L, recordForSS2.clusteredVariantAccession)

        // Ensure no action of SS2 being assigned rs2 above is recorded
        def opRecordForSS2MapWt = dbEnv.mongoTemplate.find(query(where("inactiveObjects.hashedMessage")
                .is(recordForSS2.hashedMessage).and("_id").regex("RESTORE_MAPWT_.*")), dbsnpSvoeClass)
        assertEquals(0, opRecordForSS2MapWt.size())
    }
}
