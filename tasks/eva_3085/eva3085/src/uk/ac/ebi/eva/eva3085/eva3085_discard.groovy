// This script discards SVEs and CVEs incorrectly remapped to contig CM008173.2 and assembly GCA_000247795.2
package uk.ac.ebi.eva.eva3085

import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.eva3085.EVACursor
import uk.ac.ebi.eva.eva3085.EVADatabaseEnvironment
import uk.ac.ebi.eva.eva3085.EVALoggingUtils

import static uk.ac.ebi.eva.eva3085.EVADatabaseEnvironment.*

def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, Application.class)
def scriptLogger = EVALoggingUtils.getLogger(this.class)

def assembly = "GCA_000247795.2"
def contig = "CM008173.2"

SubmittedVariantOperationEntity getDiscardSvoe(SubmittedVariantEntity sve) {
    String svoeId = "DISCARD_SS_${sve.getAccession()}_FROM_${assembly}"
    SubmittedVariantOperationEntity svoe = new SubmittedVariantOperationEntity();
    svoe.fill(EventType.DISCARDED,
              sve.getAccession(),
              "Submitted variant discarded due to being remapped to invalid contig",
              [new SubmittedVariantInactiveEntity(sve)])
    svoe.setId(svoeId)
    return svoe
}

ClusteredVariantOperationEntity getDiscardCvoe(ClusteredVariantEntity cve) {
    String cvoeId = "DISCARD_RS_${cve.getAccession()}_FROM_${assembly}"
    ClusteredVariantOperationEntity cvoe = new ClusteredVariantOperationEntity();
    cvoe.fill(EventType.DISCARDED,
              cve.getAccession(),
              "Clustered variant discarded due to being remapped to invalid contig",
              [new ClusteredVariantInactiveEntity(cve)])
    cvoe.setId(cvoeId)
    return cvoe
}


// Discard remapped submitted variants
def svesDiscarded = 0
def svesToDiscard = new EVACursor(
    where("seq").is(assembly).and("contig").is(contig).and("remappedFrom").exists(true), prodEnv.mongoTemplate, sveClass)
svesToDiscard.each{svesInBatch ->
    // Create and insert discard operations
    def discardOperations = svesInBatch.collect{sve -> getDiscardSvoe(sve)}
    prodEnv.bulkInsertIgnoreDuplicates(discardOperations, svoeClass)
    // Remove SVE documents
    def removedDocs = prodEnv.mongoTemplate.findAllAndRemove(
        query(where("_id").in(svesInBatch.collect{it.getId()})), sveClass)
    svesDiscarded += removedDocs.size()
    scriptLogger.info("${svesDiscarded} SVEs discarded thus far...")
}

// Discard clustered variants
def cvesDiscarded = 0
def cvesToDiscard = new EVACursor(where("asm").is(assembly).and("contig").is(contig), prodEnv.mongoTemplate, cveClass)
cvesToDiscard.each{cvesInBatch ->
    // Create and insert discard operations
    def discardOperations = cvesInBatch.collect{cve -> getDiscardCvoe(cve)}
    prodEnv.bulkInsertIgnoreDuplicates(discardOperations, cvoeClass)
    // Remove CVE documents
    def removedDocs = prodEnv.mongoTemplate.findAllAndRemove(
        query(where("_id").in(cvesInBatch.collect{it.getId()})), cveClass)
    cvesDiscarded += removedDocs.size()
    scriptLogger.info("${cvesDiscarded} CVEs discarded thus far...")
}

// Check that there are none left on this contig and assembly
scriptLogger.info("Remaining SVEs: ", prodEnv.mongoTemplate.count(
    query(where("seq").is(assembly).and("contig").is(contig)), sveClass))
scriptLogger.info("Remaining CVEs: ", prodEnv.mongoTemplate.count(
    query(where("asm").is(assembly).and("contig").is(contig)), cveClass))
