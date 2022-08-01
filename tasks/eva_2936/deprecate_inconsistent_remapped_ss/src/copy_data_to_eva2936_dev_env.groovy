package src

@Grab(group = 'org.apache.commons', module = 'commons-csv', version = '1.9.0')
@Grab(group='org.codehaus.groovy', module='groovy-cli-commons', version='3.0.10')

import org.apache.commons.csv.*
import org.slf4j.LoggerFactory
import src.EVADatabaseEnvironment
import src.EVADataSet
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.deprecate.Application


def cli = new CliBuilder(usage: 'copy_data_to_eva2936_dev_env <TSV folder with duplicates report> <dev properties file> <prod properties file>')
def options = cli.parse(args)
if (options.arguments().size() != 3) {
    cli.usage()
    return
}

def duplicateSSReportDir = options.arguments()[0]
def devPropertiesFile = options.arguments()[1]
def prodPropertiesFile = options.arguments()[2]

def prodEnv = EVADatabaseEnvironment.createFromSpringContext(prodPropertiesFile, Application.class)
def eva2936DevEnv = EVADatabaseEnvironment.createFromSpringContext(devPropertiesFile, Application.class)
def evaSSAccessionsLowerBound = eva2936DevEnv.springApplicationContext.getBean("accessioningMonotonicInitSs", Long.class)
def evaRSAccessionsLowerBound = eva2936DevEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class)
def logger = LoggerFactory.getLogger(this.class)
def categoriesForDbsnpDuplicateSS =
        ["In_original_assembly,Multi_position_ssid,Remapped,Same_variants_in_source",
         "In_original_assembly,Multi_allele_ssid,Remapped,Same_variants_in_source",
         "Multi_allele_ssid,Multi_position_in_source,Remapped",
         "Multi_allele_in_source,Multi_allele_ssid,Remapped",
         "Multi_position_in_source,Multi_position_ssid,Remapped",
         "In_original_assembly,Multi_position_in_source,Multi_position_ssid,Remapped",
         "Multi_allele_in_source,Multi_position_ssid,Remapped",
         "In_original_assembly,Multi_allele_in_source,Multi_position_ssid,Remapped"]
def filesWithDuplicateSS = ["dbsnp_analysis", "eva_analysis"].collect{category ->
    new FileNameByRegexFinder().getFileNames(duplicateSSReportDir + "/" + category, /.*_duplicates_annotated\.tsv/)}

// load duplicate report to Mongo - so that we don't have to grep these every time an accession needs to be looked up
filesWithDuplicateSS.each { category ->
    category.each { fileName ->
        new File(fileName).withReader { reader ->
            logger.info("Processing file ${fileName}...")
            CSVParser csv = new CSVParser(reader, CSVFormat.TDF)
            csv.findAll { csvRecord ->
                categoriesForDbsnpDuplicateSS.contains(csvRecord.get(2).trim())
            }.collect { obj -> (CSVRecord) obj }
                    .collate(1000)
                    .each { csvRecords ->
                        def docsToInsert = csvRecords.collect { csvRecord ->
                            new DuplicateSSCategory(csvRecord[0].toLong(), csvRecord[1].toString(),
                                    csvRecord[2].toString())
                        }
                        logger.info("Inserting ${docsToInsert.size()} documents...")
                        eva2936DevEnv.bulkUpsert(docsToInsert, DuplicateSSCategory.class)
                    }
        }
    }
}

// load duplicate SS to Mongo
def duplicateReports = new EVADataSet(null, eva2936DevEnv.mongoTemplate, DuplicateSSCategory.class)
duplicateReports.each {List<DuplicateSSCategory> categories ->
    def ssGroupedByAssembly = categories.groupBy {it.seq}
    ssGroupedByAssembly.each {assembly, categoriesInAssembly ->
        def ssIDs = categoriesInAssembly.collect{it.accession}
        // Groupby clause below will only look for variants that satisfy category 1 through 3 criteria here: https://docs.google.com/spreadsheets/d/1LZbUtamFtQH12SVTfCfYOjcIHWL1OOC_kYnieRLR9UA/edit#rangeid=1722186521
        // i.e., EVA remapped variants with accession collisions on existing variants from dbSNP (category 1) OR from multiple source assemblies (category 2) OR both
        def svesToInsert = prodEnv.submittedVariantAccessioningService.getAllActiveByAssemblyAndAccessionIn(assembly, ssIDs)
                .groupBy {it.getAccession()}.findAll { k, v -> return v.collect {
                                                        it.getData().getRemappedFrom() }.unique().size() > 1}
                .values().flatten().collect{new SubmittedVariantEntity(it.getAccession(), it.getHash(), it.getData(),
                it.getVersion())}
        if (svesToInsert.size() > 0) {
            def dbsnpSVEsToInsert = svesToInsert.findAll{it.getAccession() < evaSSAccessionsLowerBound}
            def evaSVEsToInsert = svesToInsert.findAll{it.getAccession() >= evaSSAccessionsLowerBound}
            logger.info("Inserting ${dbsnpSVEsToInsert.size()} to DbsnpSubmittedVariantEntity in DEV...")
            logger.info("Inserting ${evaSVEsToInsert.size()} to SubmittedVariantEntity in DEV...")
            eva2936DevEnv.bulkUpsert(dbsnpSVEsToInsert, DbsnpSubmittedVariantEntity.class)
            eva2936DevEnv.bulkUpsert(evaSVEsToInsert, SubmittedVariantEntity.class)
        }
    }
}

[SubmittedVariantEntity.class, DbsnpSubmittedVariantEntity.class].collect{collection ->
    new EVADataSet<>(null, eva2936DevEnv.mongoTemplate, collection)}.each {dataset -> dataset.each{sves ->
    sves.groupBy{it.getReferenceSequenceAccession()}.each{String assembly, svesInAssembly ->
        List<Long> rsIDs = svesInAssembly.collect{it.getClusteredVariantAccession()}
        // Insert associated RS IDs
        prodEnv.clusteredVariantAccessioningService.getAllActiveByAssemblyAndAccessionIn(assembly, rsIDs)
                .collect{new ClusteredVariantEntity(it.getAccession(), it.getHash(), it.getData(), it.getVersion())}
                .collate(1000).each{cvesToInsert ->
            def dbsnpCVEsToInsert = cvesToInsert.findAll{it.getAccession() < evaRSAccessionsLowerBound}
            def evaCVEsToInsert = cvesToInsert.findAll{it.getAccession() >= evaRSAccessionsLowerBound}
            logger.info("Inserting ${dbsnpCVEsToInsert.size()} associated RS to DbsnpClusteredVariantEntity in DEV...")
            logger.info("Inserting ${evaCVEsToInsert.size()} associated RS to ClusteredVariantEntity in DEV...")
            eva2936DevEnv.bulkUpsert(dbsnpCVEsToInsert, DbsnpClusteredVariantEntity.class)
            eva2936DevEnv.bulkUpsert(evaCVEsToInsert, ClusteredVariantEntity.class)
        }
        // Insert other SS IDs with these RS IDs
        prodEnv.submittedVariantAccessioningService.getByClusteredVariantAccessionIn(rsIDs)
                .findAll { it.getData().getReferenceSequenceAccession().equals(assembly)}
                .collect{new SubmittedVariantEntity(it.getAccession(), it.getHash(), it.getData(), it.getVersion())}
                .collate(1000).each {svesToInsert ->
            def dbsnpSVEsToInsert = svesToInsert.findAll{it.getAccession() < evaSSAccessionsLowerBound}
            def evaSVEsToInsert = svesToInsert.findAll{it.getAccession() >= evaSSAccessionsLowerBound}
            logger.info("Inserting ${dbsnpSVEsToInsert.size()} associated SS to DbsnpSubmittedVariantEntity in DEV...")
            logger.info("Inserting ${evaSVEsToInsert.size()} associated SS to SubmittedVariantEntity in DEV...")
            eva2936DevEnv.bulkUpsert(dbsnpSVEsToInsert, DbsnpSubmittedVariantEntity.class)
            eva2936DevEnv.bulkUpsert(evaSVEsToInsert, SubmittedVariantEntity.class)
        }
}}}

// Backup collections before remediation
eva2936DevEnv.mongoTemplate.getCollectionNames().each {collectionName ->
    def before_remediation_suffix = "_before_remediation"
    if (!collectionName.endsWith(before_remediation_suffix)) {
        eva2936DevEnv.mongoTemplate.getCollection(collectionName).aggregate(
                Collections.singletonList(new org.bson.Document("\$out",
                        "${collectionName}_${before_remediation_suffix}".toString()))).allowDiskUse(true).size()
    }
}
