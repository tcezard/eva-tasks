@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-clustering', version = '0.6.4')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.4')
@Grab(group = 'org.apache.commons', module = 'commons-csv', version = '1.9.0')
@Grab(group='org.codehaus.groovy', module='groovy-cli-commons', version='3.0.10')
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import org.springframework.data.mongodb.core.query.Query
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity

import static org.springframework.data.mongodb.core.query.Criteria.where
import EVADatabaseEnvironment

class LoadDuplicateDbsnpSSToTestDB {
    static void main(String[] args) {

        var cli = new CliBuilder(usage: 'LoadDuplicateDbsnpSSToTestDB <TSV folder with duplicates> <production properties file> <development properties file>')
        var options = cli.parse(args)
        if (options.arguments().size() != 3) {
            cli.usage()
            return
        }

        var duplicateSSFolder = options.arguments()[0]
        var prodPropertiesFile = options.arguments()[1]
        var devPropertiesFile = options.arguments()[2]
        EVADatabaseEnvironment prodEnv = EVADatabaseEnvironment.parseFrom(prodPropertiesFile)
        EVADatabaseEnvironment devEnv = EVADatabaseEnvironment.parseFrom(devPropertiesFile)

        def filesWithDuplicateSS = new FileNameByRegexFinder()
                .getFileNames(duplicateSSFolder, /.*_duplicates_annotated\.tsv/)

        filesWithDuplicateSS.each {fileWithDuplicateSS ->
            loadDuplicateDbsnpSSToTestDB(fileWithDuplicateSS, prodEnv, devEnv)
        }

    }

    static void loadDuplicateDbsnpSSToTestDB(String fileWithDuplicateSS, EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv) {
        println("Processing file ${fileWithDuplicateSS} to insert dbSNP SS duplicates..")
        /**
         * Parse a TSV file with missing RS that looks like this (SS ID, Source, Assembly, Category of duplicate - EVA-2840)
         208966003	GCA_002880775.3	Multi_position_in_source,Multi_position_ssid,Remapped
         */       // Categories for dbSNP duplicates - see https://www.ebi.ac.uk/panda/jira/browse/EVA-2840?focusedCommentId=392526&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-392526
        var categoriesForDbsnpDuplicateSS =
                ["In_original_assembly,Multi_allele_ssid", "In_original_assembly,Multi_position_ssid",
                 "In_original_assembly,Multi_allele_in_source,Multi_position_ssid,Remapped",
                 "In_original_assembly,Multi_position_in_source,Multi_position_ssid,Remapped"]
        new File(fileWithDuplicateSS).withReader {
            reader ->
                CSVParser csv = new CSVParser(reader, CSVFormat.TDF)
                csv.findAll  {csvRecord ->
                    categoriesForDbsnpDuplicateSS.contains(csvRecord.get(2).trim())}
                        .collect {obj -> (CSVRecord)obj}
                        .collate(1000)
                        .each { rowsInFileWithDuplicateDbsnpSS ->
                            insertDuplicateSSWithAccessionsFromProd(rowsInFileWithDuplicateDbsnpSS, prodEnv, devEnv)
                        }
        }
    }

    static void insertDuplicateSSWithAccessionsFromProd(List<CSVRecord> rowsInFileWithDuplicateDbsnpSS, EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv) {
        var accessionsByAssembly = rowsInFileWithDuplicateDbsnpSS.groupBy { row -> row.get(1).trim() }
        accessionsByAssembly.each{assembly, rows ->
            var accessionsToFind = rows.collect{row -> Long.parseLong(row.get(0))}
            Query queryToLookupDuplicateSS =
                    new Query(where("seq").in(assembly).and("accession").in(accessionsToFind))
            var dbsnpSVEsToInsert = prodEnv.mongoTemplate.find(queryToLookupDuplicateSS, DbsnpSubmittedVariantEntity.class)
            var dbsnpSVEsalreadyInserted = devEnv.mongoTemplate.find(queryToLookupDuplicateSS, DbsnpSubmittedVariantEntity.class)
            dbsnpSVEsToInsert = dbsnpSVEsToInsert - dbsnpSVEsalreadyInserted
            if (dbsnpSVEsToInsert.size() > 0) {
                devEnv.mongoTemplate.insert(dbsnpSVEsToInsert, DbsnpSubmittedVariantEntity.class)
            }
        }
    }
}
