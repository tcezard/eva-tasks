package eva3116

import groovy.cli.picocli.CliBuilder
import org.springframework.data.mongodb.core.query.Query
import uk.ac.ebi.eva.accession.core.GenericApplication

import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static org.springframework.data.mongodb.core.query.Criteria.where

// This script analyzes the root cause behind some SS entries having map-weight but the corresponding RS not having them
def cli = new CliBuilder()
cli.prodPropertiesFile(args:1, "Production properties file for accessioning", required: true)
cli.devPropertiesFile(args:1, "EVA-3097 (where initial assessment was carried out) development environment properties " +
        "file for accessioning", required: true)
def options = cli.parse(args)
if (!options) {
    cli.usage()
    System.exit(1)
}

def prodEnv = createFromSpringContext(options.prodPropertiesFile, GenericApplication.class)
def devEnv = createFromSpringContext(options.devPropertiesFile, GenericApplication.class)

def eva3116ImpactedRS = devEnv.mongoTemplate.find(new Query(), dbsnpCveClass, "cvesOfMapWtSSWithoutMapWt")
// 479
println(eva3116ImpactedRS.size())

def idsToSearch = eva3116ImpactedRS.collect{"EVA2850_MERGED_[0-9]+" + "_" + it.hashedMessage}
// Ensure that all the impacted RS IDs were merged into as part of EVA-2850 - see https://github.com/EBIvariation/eva-tasks/blob/33503fdbad9f89079f91f2fc0a530b0487c7aff2/tasks/eva_2850/fix_discordant_variants.py#L229
// These will be remedied by scripts written in EVA-2205
// 479
println(idsToSearch.collect {idToSearch ->
    prodEnv.mongoTemplate.find(query(where("_id").regex(idToSearch)), dbsnpCvoeClass)}.flatten().size())
