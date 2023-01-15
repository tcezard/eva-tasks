package eva3117

import org.slf4j.LoggerFactory
import uk.ac.ebi.eva.accession.core.GenericApplication
import uk.ac.ebi.eva.accession.core.model.ISubmittedVariant
import uk.ac.ebi.eva.groovy.commons.EVAObjectModelUtils
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where
import groovy.cli.picocli.CliBuilder

def cli = new CliBuilder()
cli.propertiesFile(args:1, "Properties file to use for remediation", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
def logger = LoggerFactory.getLogger(GenericApplication.class)

// Affected assemblies to check
// see category ss_wo_mapwt_sharing_mapwt_accessions in https://docs.google.com/spreadsheets/d/1vUQ0aZsZUJkQI_7CUYJ0_JN7Vd9oNxh7TazTOcpRuIE/edit#gid=0
def assemblies = ["GCA_000001215.2","GCA_000001635.5","GCA_000001515.4","GCA_000001735.1","GCA_000001895.3","GCA_000002035.2","GCA_000001895.4","GCA_000002035.3","GCA_000002175.2","GCA_000002195.1","GCA_000002255.2","GCA_000002285.1","GCA_000002285.2","GCA_000002305.1","GCA_000002315.1","GCA_000002315.3","GCA_000002655.1","GCA_000002265.1","GCA_000003025.4","GCA_000003025.6","GCA_000002775.1","GCA_000003055.3","GCA_000003195.1","GCA_000003055.5","GCA_000003205.4","GCA_000003205.6","GCA_000003745.2","GCA_000004515.2","GCA_000004515.3","GCA_000001635.6","GCA_000005575.1","GCA_000146605.2","GCA_000146605.3","GCA_000151805.2","GCA_000005425.2","GCA_000181335.3","GCA_000188235.2","GCA_000298735.1","GCA_000298735.2","GCA_000317375.1","GCA_000317765.1","GCA_000772875.3","GCA_001433935.1","GCA_000002315.2","GCA_000003205.1","GCA_000247795.2","GCA_000001635.4"]

// Jacketing routine to suppress errors like the following that arise due to unrepresentable variants imported from dbSNP
// ex: Cannot determine the type of the Variant with Reference: TAAGAGCCCTAGAGTTGACAGGGGTACTCATTGGTGTCCTTCACTTCAGCCTGCCCTGGGACGCATGAAGGGTCTTGCTAATGGATGTAATTCTGTTAGAATATCATTTTCTGAAATGGTTTTAGGGACAAGGGACAGATGTATTTTTCCTCTGTGTCACTATAACACATTGCAGAAGGGAGGAACTTGAGAGCCATAGAAGCAGGAAGACCCACTGTCTGCATTGGAGGAGCATTGGATGGCTCCCTTC...(TOTAL 828), Alternate: A
def getClusteredVariantHashSuppressErrors(ISubmittedVariant submittedVariant) {
    def result = null
    try {
        result = EVAObjectModelUtils.getClusteredVariantHash(submittedVariant)
    }
    catch(IllegalArgumentException illegalArgumentException) {
        return "-1"
    }
    return result
}

assemblies.each {assembly ->
def ssWithMapWtCursor = new RetryableBatchingCursor<>(where("seq").is(assembly).and(
        "mapWeight").exists(true), prodEnv.mongoTemplate, dbsnpSveClass)
def numEntriesScanned = 0
ssWithMapWtCursor.each{svesWithMapWt ->
    def svesWithoutMapWtSharingID = prodEnv.mongoTemplate.find(query(where("seq").is(assembly).and(
            "accession").in(svesWithMapWt.collect{it.accession}).and("mapWeight").exists(false)),
            dbsnpSveClass).groupBy {it.accession}
    def svesWithMapWtGroupedBySSID = svesWithMapWt.groupBy{it.accession}
    svesWithoutMapWtSharingID.each{ssID, svesWithoutMapWt ->
        def rsHashesOfSvesWithMapWt = svesWithMapWtGroupedBySSID.get(ssID).collect{getClusteredVariantHashSuppressErrors(it)}
        svesWithoutMapWt.each{sveWithoutMapWt ->
            if (rsHashesOfSvesWithMapWt.contains(getClusteredVariantHashSuppressErrors(sveWithoutMapWt))) {
                logger.error("Non map-weight SS with hash ${sveWithoutMapWt.hashedMessage} and " +
                        "accession ${sveWithoutMapWt.accession} shares ID and RS locus with a map-weighted SS")
            }
        }
    }
    numEntriesScanned += svesWithMapWt.size()
    logger.info("Scanned ${numEntriesScanned} entries in ${assembly} so far....")
}}
