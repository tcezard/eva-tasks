package src

@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-clustering', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'variation-commons-batch', version = '0.8.1')
@Grab(group = 'org.springframework.batch', module = 'spring-batch-test', version = '4.1.0.RELEASE')
@Grab(group = 'org.apache.commons', module = 'commons-csv', version = '1.9.0')
@Grab(group='org.codehaus.groovy', module='groovy-cli-commons', version='3.0.10')

import com.mongodb.WriteConcern
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import org.bson.Document
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.batch.test.MetaDataInstanceFactory
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.beans.factory.annotation.Value
import org.springframework.boot.CommandLineRunner
import org.springframework.boot.SpringApplication
import org.springframework.boot.autoconfigure.SpringBootApplication
import org.springframework.context.ConfigurableApplicationContext
import org.springframework.context.annotation.Import
import org.springframework.dao.DataIntegrityViolationException
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Query
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.accession.clustering.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.io.SSSplitWriterConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.listeners.ListenersConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.service.nonhuman.SubmittedVariantAccessioningService
import uk.ac.ebi.eva.commons.batch.io.AggregatedVcfLineMapper
import uk.ac.ebi.eva.commons.batch.io.VcfReader
import uk.ac.ebi.eva.commons.core.models.Aggregation
import uk.ac.ebi.eva.metrics.configuration.MetricConfiguration

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query

@Component
@Import(value=[SSSplitWriterConfiguration.class,
        InputParametersConfiguration.class, MetricConfiguration.class, ListenersConfiguration.class])
class LoadDuplicateRemappedSSToTestDB implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(LoadDuplicateRemappedSSToTestDB.class)

    @Autowired
    private SubmittedVariantAccessioningService submittedVariantAccessioningService

    @Autowired
    private MongoTemplate prodMongoTemplate

    private MongoTemplate devMongoTemplate

    @Value('${devenv.properties}')
    private String devEnvPropertiesFile

    @Value('${duplicate.ss.folder}')
    private String duplicateSSFolder

    @Value('${remapping.root.folder}')
    private String remappingRootFolder
    private String remappedAttributesCollection

    void run(String... args) {
        this.devMongoTemplate = EVADatabaseEnvironment.parseFrom(devEnvPropertiesFile).mongoTemplate
        this.devMongoTemplate.setWriteConcern(WriteConcern.MAJORITY)

        def filesWithDuplicateSS = new FileNameByRegexFinder()
                .getFileNames(duplicateSSFolder, /.*_duplicates_annotated\.tsv/)

        filesWithDuplicateSS.each {fileWithDuplicateSS ->
            insertDuplicateRemappedSSToTestDB(fileWithDuplicateSS)
        }

    }

    void insertDuplicateRemappedSSToTestDB(String fileWithDuplicateSS) {
        println("Processing file ${fileWithDuplicateSS} to insert dbSNP SS duplicates..")
        /**
         * Parse a TSV file with missing RS that looks like this (SS ID, Source, Assembly, Category of duplicate - EVA-2840)
         208966003	GCA_002880775.3	Multi_position_in_source,Multi_position_ssid,Remapped
         */
        // Categories for remapped SS duplicates - see categories 1 through 4 here: https://www.ebi.ac.uk/panda/jira/browse/EVA-2889?focusedCommentId=396121&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-396121
        def categoriesForDbsnpDuplicateSS =
                ["Multi_allele_ssid,Remapped,Same_variants_in_source",
                 "Multi_position_ssid,Remapped,Same_variants_in_source"]
        new File(fileWithDuplicateSS).withReader {
            reader ->
                CSVParser csv = new CSVParser(reader, CSVFormat.TDF)
                csv.findAll  {csvRecord ->
                    categoriesForDbsnpDuplicateSS.contains(csvRecord.get(2).trim())}
                        .collect {obj -> (CSVRecord)obj}
                        .collate(1000)
                        .each { rowsInFileWithDuplicateRemappedSS ->
                            insertSSRecordsToDev(rowsInFileWithDuplicateRemappedSS)
                        }
        }
    }

    void insertSSRecordsToDev(List<CSVRecord> rowsInFileWithDuplicateDbsnpSS) {
        def accessionsByAssembly = rowsInFileWithDuplicateDbsnpSS.groupBy { row -> row.get(1).trim() }
        accessionsByAssembly.each{assembly, rows ->
            def accessionsToFind = rows.collect{row -> Long.parseLong(row.get(0))}
            Query queryToLookupDuplicateRemappedSS =
                    new Query(where("seq").is(assembly).and("accession").in(accessionsToFind))
            List<DbsnpSubmittedVariantEntity> dbsnpSVEsToInsert =
                    this.prodMongoTemplate.find(queryToLookupDuplicateRemappedSS, DbsnpSubmittedVariantEntity.class)
            def dbsnpSVEsToInsertGroupedBySourceAssembly =
                    dbsnpSVEsToInsert.findAll {dbsnpSVE -> Objects.nonNull(dbsnpSVE.getRemappedFrom())}.groupBy {dbsnpSVE -> dbsnpSVE.getRemappedFrom()}
            // Remapped SS which have no corresponding entries in the original assembly in the dbsnpSVE collection
            // This probably means that these SS IDs were split during EVA-2861 and now reside in the SVE collection.
            // This shouldn't happen but we check anyway!! See categories - https://www.ebi.ac.uk/panda/jira/browse/EVA-2889?focusedCommentId=396121&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-396121
            // If we find anything in this category, they should be dealt with by EVA-2862
            def dbsnpSVEsWithNoCorrespondingEntriesInSourceAssembly =
                    getDbsnpSVEsWithNoCorrespondingEntriesInSourceAssembly(dbsnpSVEsToInsertGroupedBySourceAssembly)

            List<SubmittedVariantEntity> svesToInsert =
                    this.prodMongoTemplate.find(queryToLookupDuplicateRemappedSS, SubmittedVariantEntity.class)

            List<String> originalAssemblies = dbsnpSVEsToInsert.findAll {Objects.nonNull(it.getRemappedFrom())}.collect{it.getRemappedFrom()}
            originalAssemblies.addAll(svesToInsert.findAll {Objects.nonNull(it.getRemappedFrom())}.collect{it.getRemappedFrom()})
            originalAssemblies = originalAssemblies.unique()
            Query queryToLookupOriginalSS =
                    new Query(where("seq").in(originalAssemblies).and("accession").in(accessionsToFind))
            dbsnpSVEsToInsert.addAll(this.prodMongoTemplate.find(queryToLookupOriginalSS, DbsnpSubmittedVariantEntity.class))
            svesToInsert.addAll(this.prodMongoTemplate.find(queryToLookupOriginalSS, SubmittedVariantEntity.class))

            Map<String, List<? extends SubmittedVariantEntity>> sveAndDbsnpSVEGroupedByTaxonomyAndSeq =
                    ((dbsnpSVEsToInsert + svesToInsert) as List<? extends SubmittedVariantEntity>).findAll{Objects.nonNull(it.getRemappedFrom())}.groupBy{"" + "${it.getTaxonomyAccession()}/${it.getRemappedFrom()}"}
            insertLinesFromRemappedVCF(sveAndDbsnpSVEGroupedByTaxonomyAndSeq)

            def dbsnpSVEsAlreadyInserted = this.devMongoTemplate.find(queryToLookupDuplicateRemappedSS, DbsnpSubmittedVariantEntity.class)
            dbsnpSVEsAlreadyInserted.addAll(this.devMongoTemplate.find(queryToLookupOriginalSS, DbsnpSubmittedVariantEntity.class))
            def svesAlreadyInserted = this.devMongoTemplate.find(queryToLookupDuplicateRemappedSS, SubmittedVariantEntity.class)
            svesAlreadyInserted.addAll(this.devMongoTemplate.find(queryToLookupOriginalSS, SubmittedVariantEntity.class))

            dbsnpSVEsToInsert = (dbsnpSVEsToInsert - dbsnpSVEsAlreadyInserted) - dbsnpSVEsWithNoCorrespondingEntriesInSourceAssembly
            svesToInsert = svesToInsert - svesAlreadyInserted

            if (dbsnpSVEsToInsert.size() > 0) {
                this.devMongoTemplate.insert(dbsnpSVEsToInsert, DbsnpSubmittedVariantEntity.class)
            }
            if (svesToInsert.size() > 0) {
                this.devMongoTemplate.insert(svesToInsert, SubmittedVariantEntity.class)
            }
        }
    }

    void insertLinesFromRemappedVCF(Map<String, List<? extends SubmittedVariantEntity>> sveAndDbsnpSVEGroupedByTaxonomyAndSeq) {
        sveAndDbsnpSVEGroupedByTaxonomyAndSeq.each { taxAndSeq, sveAndDbsnpSVE ->
            def out = new StringBuffer()
            def err = new StringBuffer()
            File tempPatternFile = File.createTempFile("pattern", ".txt")
            File tempOutputFile = File.createTempFile("output", ".txt")
            File tempVCFFile = File.createTempFile("tempVCF", ".vcf")

            String assembly = taxAndSeq.split("/")[1]
            String ssToGrep = sveAndDbsnpSVE.collect { "ss" + it.getAccession() }.unique().join("\n")
            tempPatternFile.write(ssToGrep)
            def commandToGrepForSSInRemappedVCF = ["bash", "-c", '''grep -w -f ''' + tempPatternFile.getAbsolutePath() + " " + this.remappingRootFolder + "/" + taxAndSeq + "/*/*_*_remapped.vcf" + " > ${tempOutputFile.getAbsolutePath()}"]
            logger.info("Running command $commandToGrepForSSInRemappedVCF...")
            def process = commandToGrepForSSInRemappedVCF.execute()
            process.waitForProcessOutput(out, err)
            if (process.exitValue() != 0) {
                logger.error("Command $commandToGrepForSSInRemappedVCF exited with exit code ${process.exitValue()}!!")
                return
            }

            tempOutputFile.readLines().each {lineFromRemappedVCF ->
                if (!(lineFromRemappedVCF.trim().equals(""))) {
                    tempVCFFile.append(lineFromRemappedVCF.split(":")[1..-1].join(":") + ";AF=1,2,3,4,5,6;AC=1,2,3,4,5,6;\n")
                }
            }
            AggregatedVcfLineMapper lineMapper = new AggregatedVcfLineMapper("fileId", "studyId", Aggregation.BASIC, null)
            lineMapper.setIncludeIds(true)
            VcfReader vcfReader = new VcfReader(lineMapper, tempVCFFile)
            vcfReader.open(MetaDataInstanceFactory.createStepExecution().getExecutionContext())
            List<Document> variantDocumentsToInsert = new ArrayList<>()

            while(true) {
                remappedAttributesCollection = "remappedAttributes"
                def variants = vcfReader.read()
                if (Objects.isNull(variants)) break
                variants.each {variant ->
                    try {
                        Long accession = Long.parseLong(variant.getIds()[0].replace("ss", ""))
                        def chr = variant.getChromosome()
                        def start = variant.getStart()
                        def end = variant.getEnd()
                        def ref = variant.getReference()
                        def alt = variant.getAlternate()
                        String id = "${assembly}_${accession}_${chr}_${start}_${end}_${ref}_${alt}"
                        Document variantDocument = new Document("_id", id)
                        variantDocument.put("accession", accession)
                        variantDocument.put("seq", assembly)
                        variantDocument.put("ref", ref)
                        variantDocument.put("alt", alt)
                        variantDocument.put("contig", chr)
                        variantDocument.put("start", start)
                        variantDocument.put("end", end)
                        variantDocument.put("type", variant.getType())
                        variantDocument.putAll(variant.getSourceEntries()[0].getAttributes())
                        variantDocument.remove("src")
                        variantDocumentsToInsert.add(variantDocument)
                    }
                    catch(Exception ex) {
                        logger.error(ex.getMessage())
                        ex.printStackTrace()
                    }
                }
            }
            variantDocumentsToInsert = variantDocumentsToInsert.unique{doc1, doc2 -> doc1.get("_id") <=> doc2.get("_id")}
            try {
                List<String> idsToInsert = variantDocumentsToInsert.collect{doc -> doc.get("_id").toString()}
                List<String> idsExisting = this.devMongoTemplate.find(query(where("_id").in(idsToInsert)), Document.class, remappedAttributesCollection).collect{ doc -> doc.get("_id").toString()}
                variantDocumentsToInsert.removeAll {doc ->idsExisting.contains(doc.get("_id")) }
                logger.info("Total number of documents to insert ${variantDocumentsToInsert.size()}...")
                int numElementsInserted = 0
                variantDocumentsToInsert.collate(100).each {docs ->
                    String ids = docs.collect{it.get("_id")}.join(",")
                    logger.info("Inserting ids $ids...")
                    numElementsInserted += 100
                    this.devMongoTemplate.insert(docs, remappedAttributesCollection)
                    logger.info("Inserted $numElementsInserted so far...")
                }
            }
            catch (DataIntegrityViolationException dao) {
                logger.error("Error due to ${dao.getCause().getMessage()}...")
            }
        }
    }

    List<DbsnpSubmittedVariantEntity> getDbsnpSVEsWithNoCorrespondingEntriesInSourceAssembly (
            Map<String, List<DbsnpSubmittedVariantEntity>> dbsnpSVEsToInsertGroupedBySourceAssembly) {
        return dbsnpSVEsToInsertGroupedBySourceAssembly.collect  { sourceAssembly, dbsnpSVEs ->
            def accessionsFromRemappedAssembly = dbsnpSVEs.collect{dbsnpSVE -> dbsnpSVE.getAccession()}.unique()
            Query queryToLookupSourceAssemblySS =
                    new Query(where("seq").is(sourceAssembly).and("accession").in(accessionsFromRemappedAssembly))
            List<Long> accessionsInSourceAssembly =
                    this.prodMongoTemplate.find(queryToLookupSourceAssemblySS, DbsnpSubmittedVariantEntity.class)
                            .collect{dbsnpSVEInSourceAssembly -> dbsnpSVEInSourceAssembly.getAccession()}
                            .unique()
            def accessionsNotInSourceAssembly = accessionsFromRemappedAssembly - accessionsInSourceAssembly
            accessionsNotInSourceAssembly.each {accessionNotInSourceAssembly ->
                println("ERROR: No entries could be found in the source assembly " +
                        "for the remapped variants with accession ${accessionNotInSourceAssembly}!! " +
                        "This should NOT happen!!")
            }
            return accessionsNotInSourceAssembly
        }.flatten() as List<DbsnpSubmittedVariantEntity>
    }
}

@SpringBootApplication
class LoadDuplicateRemappedSSToTestDBApp  {

    static void main(String[] args) {
        ConfigurableApplicationContext context = SpringApplication.run(LoadDuplicateRemappedSSToTestDBApp.class, args)
        System.exit(SpringApplication.exit(context))
    }

}