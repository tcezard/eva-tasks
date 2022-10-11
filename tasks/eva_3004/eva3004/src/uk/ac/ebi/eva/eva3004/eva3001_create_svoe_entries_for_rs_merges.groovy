// This script creates corresponding SS merge entries for RS merges that were carried out in EVA-2850
// See here: https://ebi-eva.slack.com/archives/C01BKSCFBEG/p1664983447926359

package uk.ac.ebi.eva.eva3004

import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.eva3004.EVADatabaseEnvironment
import uk.ac.ebi.eva.eva3004.EVALoggingUtils

import static uk.ac.ebi.eva.eva3004.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning",  required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class)
def scriptLogger = EVALoggingUtils.getLogger(this.class)

def allEVA2850MergeOps = [cvoeClass, dbsnpCvoeClass].collect{prodEnv.mongoTemplate.find(query(where("_id")
        .regex("EVA2850_MERGED.*")), it)}.flatten()
allEVA2850MergeOps.collate(1000).each {ops ->
    ops.groupBy{it.inactiveObjects[0].assemblyAccession}.each {assembly, opsInBatch ->
        def rsIDsToLookUp = opsInBatch.collect { it.mergedInto }
        def opsGroupedByRSHash = opsInBatch.groupBy {it.inactiveObjects[0].hashedMessage}
        def assignedSSGroupedByRSHash = [sveClass, dbsnpSveClass].collect{prodEnv.mongoTemplate.find(query(
                where("rs").in(rsIDsToLookUp).and("seq").is(assembly)), it)}.flatten().groupBy{
            EVADatabaseEnvironment.toClusteredVariantEntity(it).hashedMessage}
        (opsGroupedByRSHash.keySet() - assignedSSGroupedByRSHash.keySet()).each {
            scriptLogger.warn("SS with RS hash ${it} and IDs ${opsGroupedByRSHash[it].mergedInto} could not be located...")
        }
        def svoeOps = assignedSSGroupedByRSHash.collect{hash, svesWithHash ->
            def clusteredOpsWithHash = opsGroupedByRSHash[hash]
            clusteredOpsWithHash.collect{clusteredOp -> svesWithHash.collect { sveWithHash ->
                sveWithHash.clusteredVariantAccession = clusteredOp.accession
                def svoeOp = new SubmittedVariantOperationEntity()
                svoeOp.fill(EventType.UPDATED, sveWithHash.accession,
                        "Original rs${clusteredOp.accession} was merged into rs${clusteredOp.mergedInto}.",
                        Arrays.asList(new SubmittedVariantInactiveEntity(sveWithHash)))
                svoeOp.setId("EVA2850_MERGED_${sveWithHash.accession}_${sveWithHash.hashedMessage}")
                return svoeOp
            }}
        }.flatten()
        prodEnv.bulkInsertIgnoreDuplicates(svoeOps.findAll{it.accession < 5e9}, dbsnpSvoeClass)
        prodEnv.bulkInsertIgnoreDuplicates(svoeOps.findAll{it.accession >= 5e9}, svoeClass)
    }
}

/**
 * Validation of the results from the above script
 */
def svoeRSAndAssembly = [dbsnpSvoeClass, svoeClass].collect{prodEnv.mongoTemplate.find(query(where("_id").regex("EVA2850_MERGED")),
        it)}.flatten().collect{svoe -> return "" + "${svoe.inactiveObjects[0].referenceSequenceAccession}_${svoe.inactiveObjects[0].clusteredVariantAccession}"}.toSet()
def cvoeRSAndAssembly = allEVA2850MergeOps.collect{cvoe -> return ""+ "${cvoe.inactiveObjects[0].assemblyAccession}_${cvoe.accession}"}.toSet()

// Ensure following does not print anything
cvoeRSAndAssembly.collate(1000).each{cvoeRSAndAssemblyInBatch ->
    (cvoeRSAndAssemblyInBatch - svoeRSAndAssembly).each{println(it)}
}
