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
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

@RunWith(SpringRunner.class)
@TestPropertySource("classpath:application-test.properties")
public class DeprecateMapWtSSTest {
    @Autowired
    private ApplicationContext applicationContext

    private boolean dbEnvSetupDone = false

    private static final String ASSEMBLY = "GCA_000000001.1"

    private static final int TAXONOMY = 60711

    private static EVADatabaseEnvironment dbEnv

    def rsLocus1, rsLocus2, ss1, ss2
    def rs1WithLocus1, rs1WithLocus2


    private static void cleanup() {
        dbEnv.mongoClient.dropDatabase(dbEnv.mongoTemplate.db.name)
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

    private void setupDB() {
        rsLocus1  = ["start": 100L, "ref": "A", "alt": "C"]
        rsLocus2  = ["start": 101L, "ref": "T", "alt": "G"]
        ss1 = createSS(1, 1, rsLocus1, 2)
        ss2 = createSS(2, 1, rsLocus2, 2)
        (rs1WithLocus1, rs1WithLocus2) =
                [ss1, ss2].collect {EVAObjectModelUtils.toClusteredVariantEntity(it) }
        dbEnv.bulkInsertIgnoreDuplicates([ss1, ss2], dbsnpSveClass)
        dbEnv.bulkInsertIgnoreDuplicates([rs1WithLocus1, rs1WithLocus2], dbsnpCveClass)
    }

    @Test
    void deprecateSS() {
        setupDB()
        DeprecateMapWtSS.deprecateSS(dbEnv, Arrays.asList(ss1, ss2))
        assertEquals(0, dbEnv.mongoTemplate.findAll(dbsnpSveClass).size())
        assertEquals(0, dbEnv.mongoTemplate.findAll(dbsnpCveClass).size())
    }
}