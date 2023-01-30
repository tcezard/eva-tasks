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

import uk.ac.ebi.eva.accession.core.model.SubmittedVariant
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.CommonTestUtils
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.EVAObjectModelUtils

import static org.junit.Assert.assertEquals
import static org.junit.Assert.assertTrue
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

@RunWith(SpringRunner.class)
@TestPropertySource("classpath:application-test.properties")
class UndoMergesIntoMultimapsTest {
    @Autowired
    private ApplicationContext applicationContext

    private boolean dbEnvSetupDone = false

    private static final String ASSEMBLY = "GCA_000000001.1"

    private static final int TAXONOMY = 60711

    private static EVADatabaseEnvironment dbEnv

    private static List<SubmittedVariantOperationEntity> svoeMerges

    def rsLocus1, rsLocus2, ss1, ss2, ss2WithRS4, ss2WithRS2, ss2WithRS3
    def rs1WithLocus1, rs1WithLocus2, rs4WithLocus2, rs2WithLocus2, rs3WithLocus2
    def svoeOp1, svoeOp2, svoeOp3
    def cvoeOp1, cvoeOp2, cvoeOp3

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
            dbEnv = createFromSpringContext(dynamicPropertyFilePath, Application.class,
                    ["parameters.assemblyAccession": ASSEMBLY])
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

    private static def mergeCvoeOp(ClusteredVariantEntity impactedRS, Long destRSID) {
        def op = new ClusteredVariantOperationEntity()
        op.fill(EventType.MERGED, impactedRS.accession, destRSID,
                "Identical clustered variant received multiple RS identifiers.",
                Arrays.asList(new ClusteredVariantInactiveEntity(impactedRS)))
        return op
    }

    private void setupScenario1() {
        rsLocus1  = ["start": 100L, "ref": "A", "alt": "C"]
        rsLocus2  = ["start": 101L, "ref": "T", "alt": "G"]
        ss1 = createSS(1, 1, rsLocus1, 2)
        ss2 = createSS(2, 1, rsLocus2, 2)
        ss2WithRS4 = createSS(2, 4, rsLocus2)
        ss2WithRS2 = createSS(2, 2, rsLocus2)
        ss2WithRS3 = createSS(2, 3, rsLocus2, 3)
        (rs1WithLocus1, rs1WithLocus2, rs4WithLocus2, rs2WithLocus2, rs3WithLocus2) =
        [ss1, ss2, ss2WithRS4, ss2WithRS2, ss2WithRS3].collect {EVAObjectModelUtils.toClusteredVariantEntity(it) }

        svoeOp1 = mergeSvoeOp(ss2WithRS4, 2)
        svoeOp2 = mergeSvoeOp(ss2WithRS2, 3)
        svoeOp3 = mergeSvoeOp(ss2WithRS3, 1)

        cvoeOp1 = mergeCvoeOp(rs4WithLocus2, rs2WithLocus2.accession)
        cvoeOp2 = mergeCvoeOp(rs2WithLocus2, rs3WithLocus2.accession)
        cvoeOp3 = mergeCvoeOp(rs3WithLocus2, rs1WithLocus2.accession)

        dbEnv.bulkInsertIgnoreDuplicates([ss1, ss2], dbsnpSveClass)
        dbEnv.bulkInsertIgnoreDuplicates([rs1WithLocus1, rs1WithLocus2], dbsnpCveClass)
        dbEnv.bulkInsertIgnoreDuplicates([svoeOp1, svoeOp2, svoeOp3], dbsnpSvoeClass)
        dbEnv.bulkInsertIgnoreDuplicates([cvoeOp1, cvoeOp2, cvoeOp3], dbsnpCvoeClass)
    }

    private void setupScenario2() {
        rsLocus1  = ["start": 100L, "ref": "A", "alt": "C"]
        rsLocus2  = ["start": 101L, "ref": "T", "alt": "G"]
        ss1 = createSS(1, 1, rsLocus1, 2)
        ss2 = createSS(2, 1, rsLocus2, 2)
        ss2WithRS4 = createSS(2, 4, rsLocus2, 3)
        ss2WithRS2 = createSS(2, 2, rsLocus2)
        ss2WithRS3 = createSS(2, 3, rsLocus2, 3)
        (rs1WithLocus1, rs1WithLocus2, rs4WithLocus2, rs2WithLocus2, rs3WithLocus2) =
                [ss1, ss2, ss2WithRS4, ss2WithRS2, ss2WithRS3].collect {EVAObjectModelUtils.toClusteredVariantEntity(it) }

        svoeOp1 = mergeSvoeOp(ss2WithRS4, 2)
        svoeOp2 = mergeSvoeOp(ss2WithRS2, 3)
        svoeOp3 = mergeSvoeOp(ss2WithRS3, 1)

        cvoeOp1 = mergeCvoeOp(rs4WithLocus2, rs2WithLocus2.accession)
        cvoeOp2 = mergeCvoeOp(rs2WithLocus2, rs3WithLocus2.accession)
        cvoeOp3 = mergeCvoeOp(rs3WithLocus2, rs1WithLocus2.accession)

        dbEnv.bulkInsertIgnoreDuplicates([ss1, ss2], dbsnpSveClass)
        dbEnv.bulkInsertIgnoreDuplicates([rs1WithLocus1, rs1WithLocus2], dbsnpCveClass)
        dbEnv.bulkInsertIgnoreDuplicates([svoeOp1, svoeOp2, svoeOp3], dbsnpSvoeClass)
        dbEnv.bulkInsertIgnoreDuplicates([cvoeOp1, cvoeOp2, cvoeOp3], dbsnpCvoeClass)
    }

    @Test
    // See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=1814908308
    void testMergeIntoMultimapScenario1() {
        setupScenario1()
        new UndoMergesIntoMultiMaps(dbEnv, ASSEMBLY).undoMergesIntoMultiMapRS()

        // Ensure that there is only one record in dbsnpSVE i.e., ss2 updated with rs2 and no mapweight
        def dbsnpSVERecords = dbEnv.mongoTemplate.findAll(dbsnpSveClass)
        def updatedRecordForSS2 = dbsnpSVERecords[0]
        assertEquals(1, dbsnpSVERecords.size())
        assertTrue(Objects.isNull(updatedRecordForSS2.mapWeight))
        assertEquals(2L, updatedRecordForSS2.clusteredVariantAccession)

        // Ensure that there is only one record in dbsnpCVE i.e., rs2 with rsLocus2
        def dbsnpCVERecords = dbEnv.mongoTemplate.findAll(dbsnpCveClass)
        assertEquals(1, dbsnpCVERecords.size())
        assertEquals(2L, dbsnpCVERecords[0].accession)
        assertEquals(rs2WithLocus2.hashedMessage, dbsnpCVERecords[0].hashedMessage)

        // Ensure the action of SS2 being assigned rs2 above is recorded
        def svoeOpRecordForSS2Update = dbEnv.mongoTemplate.find(query(where("_id").regex(
                "UNDO_MULTIMAP_RS_ASSIGN_.*")), dbsnpSvoeClass)[0]
        assertEquals("UNDO_MULTIMAP_RS_ASSIGN_" + updatedRecordForSS2.hashedMessage, svoeOpRecordForSS2Update.id)
        assertEquals(1L, svoeOpRecordForSS2Update.inactiveObjects[0].clusteredVariantAccession)
        assertTrue(Objects.nonNull(svoeOpRecordForSS2Update.inactiveObjects[0].mapWeight))

        def cvoeOpRecordForSS2Update = dbEnv.mongoTemplate.find(query(where("_id").regex(
                "UNDO_MERGE_INTO_MULTIMAP_.*")), dbsnpCvoeClass)[0]
        assertEquals("UNDO_MERGE_INTO_MULTIMAP_" + rs1WithLocus2.hashedMessage, cvoeOpRecordForSS2Update.id)
        assertEquals(EventType.UNDO_MERGE, cvoeOpRecordForSS2Update.eventType)
        assertEquals(1L, cvoeOpRecordForSS2Update.inactiveObjects[0].accession)
        assertTrue(Objects.nonNull(cvoeOpRecordForSS2Update.inactiveObjects[0].mapWeight))

        // Ensure the action of SS1 being deprecated is recorded
        def opRecordForSS1Deprecation = dbEnv.mongoTemplate.find(query(where("_id")
                .regex("SS_DEPRECATED_EVA2205_.*")), dbsnpSvoeClass)[0]
        assertEquals("SS_DEPRECATED_EVA2205_" + ss1.hashedMessage, opRecordForSS1Deprecation.id)
        assertEquals(1L, opRecordForSS1Deprecation.inactiveObjects[0].clusteredVariantAccession)
        assertTrue(Objects.nonNull(opRecordForSS1Deprecation.inactiveObjects[0].mapWeight))

        def opRecordForRS1Locus1Deprecation = dbEnv.mongoTemplate.find(query(where("_id")
                .regex("RS_DEPRECATED_EVA2205_" + rs1WithLocus1.hashedMessage)), dbsnpCvoeClass)[0]
        assertEquals(rs1WithLocus1.hashedMessage, opRecordForRS1Locus1Deprecation.inactiveObjects[0].hashedMessage)
        assertEquals(1L, opRecordForRS1Locus1Deprecation.inactiveObjects[0].accession)
        assertTrue(Objects.nonNull(opRecordForRS1Locus1Deprecation.inactiveObjects[0].mapWeight))

        def opRecordForRS1Locus2Deprecation = dbEnv.mongoTemplate.find(query(where("_id")
                .regex("RS_DEPRECATED_EVA2205_" + rs1WithLocus2.hashedMessage)), dbsnpCvoeClass)[0]
        assertEquals(rs1WithLocus2.hashedMessage, opRecordForRS1Locus2Deprecation.inactiveObjects[0].hashedMessage)
        assertEquals(1L, opRecordForRS1Locus2Deprecation.inactiveObjects[0].accession)
        assertTrue(Objects.nonNull(opRecordForRS1Locus2Deprecation.inactiveObjects[0].mapWeight))
    }

    @Test
    // See https://docs.google.com/spreadsheets/d/1kRKHR4zrq-nxH_Sg82TwXZbxdgjB3K-q4YMyJva8yMg/edit#rangeid=1471530051
    void testMergeIntoMultimapScenario2() {
        setupScenario2()
        new UndoMergesIntoMultiMaps(dbEnv, ASSEMBLY).undoMergesIntoMultiMapRS()

        // Ensure that there are no records in dbsnpSVE and dbsnpCVE since they are all deprecated
        assertEquals(0, dbEnv.mongoTemplate.findAll(dbsnpSveClass).size())
        assertEquals(0, dbEnv.mongoTemplate.findAll(dbsnpCveClass).size())

        // Ensure the action of SS1 and SS2 being deprecated is recorded
        def opRecordForSS1Deprecation = dbEnv.mongoTemplate.find(query(where("_id")
                .regex("SS_DEPRECATED_EVA2205_" + ss1.hashedMessage)), dbsnpSvoeClass)[0]
        assertEquals(1L, opRecordForSS1Deprecation.inactiveObjects[0].clusteredVariantAccession)
        assertTrue(Objects.nonNull(opRecordForSS1Deprecation.inactiveObjects[0].mapWeight))

        def opRecordForSS2Deprecation = dbEnv.mongoTemplate.find(query(where("_id")
                .regex("SS_DEPRECATED_EVA2205_" + ss2.hashedMessage)), dbsnpSvoeClass)[0]
        assertEquals("SS_DEPRECATED_EVA2205_" + ss2.hashedMessage, opRecordForSS2Deprecation.id)
        assertEquals(1L, opRecordForSS2Deprecation.inactiveObjects[0].clusteredVariantAccession)
        assertTrue(Objects.nonNull(opRecordForSS2Deprecation.inactiveObjects[0].mapWeight))

        def opRecordForRS1Locus1Deprecation = dbEnv.mongoTemplate.find(query(where("_id")
                .regex("RS_DEPRECATED_EVA2205_" + rs1WithLocus1.hashedMessage)), dbsnpCvoeClass)[0]
        assertEquals(rs1WithLocus1.hashedMessage, opRecordForRS1Locus1Deprecation.inactiveObjects[0].hashedMessage)
        assertEquals(1L, opRecordForRS1Locus1Deprecation.inactiveObjects[0].accession)
        assertTrue(Objects.nonNull(opRecordForRS1Locus1Deprecation.inactiveObjects[0].mapWeight))

        def opRecordForRS1Locus2Deprecation = dbEnv.mongoTemplate.find(query(where("_id")
                .regex("RS_DEPRECATED_EVA2205_" + rs1WithLocus2.hashedMessage)), dbsnpCvoeClass)[0]
        assertEquals(rs1WithLocus2.hashedMessage, opRecordForRS1Locus2Deprecation.inactiveObjects[0].hashedMessage)
        assertEquals(1L, opRecordForRS1Locus2Deprecation.inactiveObjects[0].accession)
        assertTrue(Objects.nonNull(opRecordForRS1Locus2Deprecation.inactiveObjects[0].mapWeight))
    }
}
