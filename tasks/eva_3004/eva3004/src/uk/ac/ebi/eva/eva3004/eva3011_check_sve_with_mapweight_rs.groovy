/***
 * This script checks if SVE which are assigned to map-weighted RS were explicitly created by clustering
 * or were just a result of EVA-2862 (split SS due to duplicates) and EVA-2905 (propagate split SS from EVA-2862)
 */

package uk.ac.ebi.eva.eva3004

import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.eva3004.EVACursor
import uk.ac.ebi.eva.accession.deprecate.Application

import static uk.ac.ebi.eva.eva3004.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, longOpt: "prod-props-file", "Production properties file for accessioning",  required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}
def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class)

// All the assemblies in DbsnpClusteredVariantEntity collection
def allAssemblies = ["GCA_000001215.2","GCA_000001215.4","GCA_000001515.4","GCA_000001515.5","GCA_000001545.3","GCA_000001635.4","GCA_000001635.5","GCA_000001635.6","GCA_000001635.9","GCA_000001735.1","GCA_000001895.3","GCA_000001895.4","GCA_000002035.2","GCA_000002035.3","GCA_000002035.4","GCA_000002175.2","GCA_000002195.1","GCA_000002255.2","GCA_000002265.1","GCA_000002285.1","GCA_000002285.2","GCA_000002305.1","GCA_000002315.1","GCA_000002315.2","GCA_000002315.3","GCA_000002315.5","GCA_000002655.1","GCA_000002775.1","GCA_000002775.3","GCA_000003025.4","GCA_000003025.6","GCA_000003055.3","GCA_000003055.5","GCA_000003195.1","GCA_000003195.3","GCA_000003205.1","GCA_000003205.4","GCA_000003205.6","GCA_000003625.1","GCA_000003745.2","GCA_000004515.2","GCA_000004515.3","GCA_000004515.4","GCA_000004665.1","GCA_000005005.5","GCA_000005005.6","GCA_000005425.2","GCA_000005575.1","GCA_000146605.2","GCA_000146605.3","GCA_000146605.4","GCA_000146795.3","GCA_000148765.2","GCA_000151805.2","GCA_000151905.1","GCA_000151905.3","GCA_000181335.3","GCA_000181335.4","GCA_000188115.2","GCA_000188115.3","GCA_000188235.2","GCA_000219495.1","GCA_000219495.2","GCA_000224145.1","GCA_000224145.2","GCA_000233375.4","GCA_000247795.2","GCA_000247815.2","GCA_000298735.1","GCA_000298735.2","GCA_000309985.1","GCA_000317375.1","GCA_000317765.1","GCA_000331145.1","GCA_000364345.1","GCA_000409795.2","GCA_000442705.1","GCA_000512255.2","GCA_000686985.1","GCA_000695525.1","GCA_000710875.1","GCA_000751015.1","GCA_000772875.3","GCA_000987745.1","GCA_001433935.1","GCA_001465895.2","GCA_001522545.1","GCA_001522545.2","GCA_001577835.1","GCA_001625215.1","GCA_001704415.1","GCA_001858045.3","GCA_002114115.1","GCA_002263795.2","GCA_002754865.1","GCA_002863925.1","GCA_002880775.3","GCA_003254395.2","GCA_003339765.3","GCA_003957565.2","GCA_011100615.1","GCA_014441545.1","GCA_015227675.2","GCA_016772045.1","GCA_902167145.1"]
allAssemblies.each{assembly ->
    println("Processing assembly ${assembly}...")
    def dbsnpCvesWithMapWtSet = new EVACursor(where("asm").is(assembly).and("mapWeight").exists(true),
            prodEnv.mongoTemplate, dbsnpCveClass)
    dbsnpCvesWithMapWtSet.each{dbsnpCvesInBatch ->
        def rsIDsToLookup = dbsnpCvesInBatch.collect{it.accession}.toSet()
        def svesWithMapWtRS = prodEnv.mongoTemplate.find(query(
                where("seq").is(assembly).and("mapWeight").exists(false).and("rs").in(rsIDsToLookup)), sveClass)
        def ssIDsWithMapWtRS = svesWithMapWtRS.collect{it.accession}
        // We don't look for specific assembly here because SS split in EVA-2862 was propagated to remapped assemblies in EVA-2905
        // Therefore the split SS is likely to be recorded in the operations for the source assembly and NOT the remapped assembly
        def svesSplitFromDbsnpSS = prodEnv.mongoTemplate.find(query(where("_id")
                .regex("SS_SPLIT_.*").and("splitInto").in(ssIDsWithMapWtRS)), dbsnpSvoeClass)
        (ssIDsWithMapWtRS - svesSplitFromDbsnpSS.collect{it.splitInto}).each {
            println("SS ID ${it} likely clustered with a map-weighted RS...")
        }
    }
}
