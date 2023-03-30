package eva3215

import groovy.cli.picocli.CliBuilder
import uk.ac.ebi.eva.accession.deprecate.Application
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

// This script deprecates orphaned RS by type
// i.e., RS with a given RS locus not associated to any SS with that RS locus
def cli = new CliBuilder()
cli.propertiesFile(args: 1, "Production properties file to use for analysis", required: true)
cli.assemblyToAnalyze(args: 1, "Assembly to analyze", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

// Since this script might be invoked multiple times across assemblies
// limit the number of connections per invocation to just 1 to avoid contention for the Spring job repository
def dbEnv = createFromSpringContext(options.propertiesFile, Application.class,
        ["spring.datasource.hikari.maximumPoolSize": 1, "parameters.assemblyAccession": options.assemblyToAnalyze])
def numSSEntriesScanned = 0
// Check SS duplicates in accession
[dbsnpSveClass, sveClass].each{collectionClass ->
    new RetryableBatchingCursor<>(where("seq").is(options.assemblyToAnalyze), dbEnv.mongoTemplate, collectionClass).each {sves ->
        def accessions = sves.collect{it.accession}.toSet()
        dbEnv.mongoTemplate.find(query(where("seq").is(options.assemblyToAnalyze).and("accession").in(accessions)),
                collectionClass).groupBy{it.accession}.findAll{it.value.size() > 1}.each {accession, records ->
            println("ERROR: Encountered SS duplicates with accession ${accession} in ${options.assemblyToAnalyze}!!")
        }
        numSSEntriesScanned += sves.size()
        println("${numSSEntriesScanned} SS entries scanned so far in ${options.assemblyToAnalyze}...")
}}
// Check RS duplicates in accession
def numRSEntriesScanned = 0
[dbsnpCveClass, cveClass].each{collectionClass ->
    new RetryableBatchingCursor<>(where("asm").is(options.assemblyToAnalyze), dbEnv.mongoTemplate, collectionClass).each {cves ->
        def accessions = cves.collect{it.accession}.toSet()
        dbEnv.mongoTemplate.find(query(where("asm").is(options.assemblyToAnalyze).and("accession").in(accessions)),
                collectionClass).groupBy{it.accession}.findAll{it.value.size() > 1}.each {accession, records ->
            println("ERROR: Encountered RS duplicates with accession ${accession} in ${options.assemblyToAnalyze}!!")
        }
        numRSEntriesScanned += cves.size()
        println("${numRSEntriesScanned} RS entries scanned so far in ${options.assemblyToAnalyze}...")
}}

// Check SS duplicates in hash
numSSEntriesScanned = 0
[dbsnpSveClass, sveClass].each{collectionClass ->
    new RetryableBatchingCursor<>(where("seq").is(options.assemblyToAnalyze), dbEnv.mongoTemplate, collectionClass).each {sves ->
        def otherCollection = (collectionClass == dbsnpSveClass)? sveClass: dbsnpSveClass
        def hashes = sves.collect{it.hashedMessage}.toSet()
        dbEnv.mongoTemplate.find(query(where("seq").is(options.assemblyToAnalyze).and("_id").in(hashes)),
                otherCollection).each {sve ->
            println("ERROR: Encountered SS hash ${sve.hashedMessage} from ${collectionClass.simpleName} in ${otherCollection.simpleName} in ${options.assemblyToAnalyze}!!")
        }
        numSSEntriesScanned += sves.size()
        println("${numSSEntriesScanned} SS entries scanned so far in ${options.assemblyToAnalyze}...")
}}

// Check SS duplicates in hash
numRSEntriesScanned = 0
[dbsnpCveClass, cveClass].each{collectionClass ->
    new RetryableBatchingCursor<>(where("asm").is(options.assemblyToAnalyze), dbEnv.mongoTemplate, collectionClass).each {cves ->
        def otherCollection = (collectionClass == dbsnpCveClass)? cveClass: dbsnpCveClass
        def hashes = cves.collect{it.hashedMessage}.toSet()
        dbEnv.mongoTemplate.find(query(where("asm").is(options.assemblyToAnalyze).and("_id").in(hashes)),
                otherCollection).each {sve ->
            println("ERROR: Encountered RS hash ${cve.hashedMessage} from ${collectionClass.simpleName} in ${otherCollection.simpleName} in ${options.assemblyToAnalyze}!!")
        }
        numRSEntriesScanned += cves.size()
        println("${numRSEntriesScanned} RS entries scanned so far in ${options.assemblyToAnalyze}...")
}}
