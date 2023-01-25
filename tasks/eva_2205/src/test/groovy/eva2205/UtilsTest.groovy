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
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.CommonTestUtils
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment

import static eva2205.eva2205_utils.*
import static org.junit.Assert.assertEquals
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

@RunWith(SpringRunner.class)
@TestPropertySource("classpath:application-test.properties")
class UtilsTest {
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

    void setupSvoeOps() {
        svoeMerges = (1..4).collect{new SubmittedVariantOperationEntity()}
        svoeMerges[0].fill(EventType.UPDATED, 1, "Original rs1 was merged into rs2.",
                Arrays.asList(new SubmittedVariantInactiveEntity(new SubmittedVariantEntity(1L, "hash1", ASSEMBLY, TAXONOMY,
                "project1", "chr1", 100L, "A", "C", 1L, true, true, true, true, 1))))
        svoeMerges[1].fill(EventType.UPDATED, 1, "Original rs3 was merged into rs4.",
                Arrays.asList(new SubmittedVariantInactiveEntity(new SubmittedVariantEntity(1L, "hash1", ASSEMBLY, TAXONOMY,
                        "project1", "chr1", 100L, "A", "C", 3L, true, true, true, true, 1))))
        svoeMerges[2].fill(EventType.UPDATED, 1, "Original rs2 was merged into rs3.",
                Arrays.asList(new SubmittedVariantInactiveEntity(new SubmittedVariantEntity(1L, "hash1", ASSEMBLY, TAXONOMY,
                        "project1", "chr1", 100L, "A", "C", 2L, true, true, true, true, 1, 3))))
        svoeMerges[3].fill(EventType.UPDATED, 1, "Original rs4 was merged into rs5.",
                Arrays.asList(new SubmittedVariantInactiveEntity(new SubmittedVariantEntity(1L, "hash1", ASSEMBLY, TAXONOMY,
                        "project1", "chr1", 100L, "A", "C", 4L, true, true, true, true, 1, 3))))
        (0..3).each{svoeMerges[it].setId("" + it)}
    }

    @Test
    void testGetSSHistoryInvolvedInRSMerges() {
        setupSvoeOps()
        dbEnv.bulkInsertIgnoreDuplicates(svoeMerges, dbsnpSvoeClass)
        // Given just one of the SS merge operations, check if the entire merge history is retrieved
        def ssHistory = getSSHistoryInvolvedInRSMerges(dbEnv, ASSEMBLY, Arrays.asList(svoeMerges[2]))
        assertEquals(4, ssHistory.size())
    }

    @Test
    void testGetChronologicalMergeChain() {
        setupSvoeOps()
        dbEnv.bulkInsertIgnoreDuplicates(svoeMerges, dbsnpSvoeClass)
        def ssHistory = getSSHistoryInvolvedInRSMerges(dbEnv, ASSEMBLY, Arrays.asList(svoeMerges[2]))
        def mergeChain = getChronologicalMergeChain(ssHistory)
        assertEquals(4, mergeChain.size())
        assertEquals(["0", "2", "1", "3"], mergeChain.collect{it.id})
    }

    @Test
    void testMostRecentNonMapWtSSRecord() {
        setupSvoeOps()
        dbEnv.bulkInsertIgnoreDuplicates(svoeMerges, dbsnpSvoeClass)
        def ssHistory = getSSHistoryInvolvedInRSMerges(dbEnv, ASSEMBLY, Arrays.asList(svoeMerges[2]))
        assertEquals("1", mostRecentNonMapWtSSRecord(getChronologicalMergeChain(ssHistory)).id)
    }

    @Test
    void testMostRecentMapWtSSRecord() {
        setupSvoeOps()
        dbEnv.bulkInsertIgnoreDuplicates(svoeMerges, dbsnpSvoeClass)
        def ssHistory = getSSHistoryInvolvedInRSMerges(dbEnv, ASSEMBLY, Arrays.asList(svoeMerges[2]))
        assertEquals("3", mostRecentMapWtSSRecord(getChronologicalMergeChain(ssHistory)).id)
    }
}
