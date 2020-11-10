import bz2
import click
import json
import traceback

from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
from pymongo import MongoClient


logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)


updated_rs_that_should_not_exist_in_EVA = []
non_updated_rs_that_should_exist_in_EVA = []
min_rs_records_to_check = 1000


def get_mongo_query_result(rs_to_find, mongo_connection_handle):
    results = []
    dbsnpCVEHandle = mongo_connection_handle["eva_accession_human_sharded"]["dbsnpClusteredVariantEntity"]
    dbsnpCVOEHandle = mongo_connection_handle["eva_accession_human_sharded"]["dbsnpClusteredVariantOperationEntity"]
    for result in dbsnpCVEHandle.find({"accession": {"$in": rs_to_find}}):
        results.append(result)
    for result in dbsnpCVOEHandle.find({"accession": {"$in": rs_to_find}}):
        results.append(result)
    return [result["accession"] for result in results]


def ensure_existing_rs_in_eva_human_accession_db(mongo_connection_handle: MongoClient, eva_production_human_dbsnp_build,
                                                 rs_id=None):
    global non_updated_rs_that_should_exist_in_EVA, min_rs_records_to_check

    num_rs_records_to_check = min_rs_records_to_check
    if rs_id:
        non_updated_rs_that_should_exist_in_EVA.append(rs_id)
    else:
        num_rs_records_to_check = len(non_updated_rs_that_should_exist_in_EVA)
    if len(non_updated_rs_that_should_exist_in_EVA) == num_rs_records_to_check or \
            (rs_id is None and num_rs_records_to_check > 0):
        results = get_mongo_query_result(non_updated_rs_that_should_exist_in_EVA, mongo_connection_handle)
        if len(results) < num_rs_records_to_check:
            logger.error("Could not find RS IDs {0} in EVA production even though they were marked in JSON "
                         "as released on or before dbSNP build {1} currently in EVA production!!"
                         .format(set(non_updated_rs_that_should_exist_in_EVA) - set(results),
                                 eva_production_human_dbsnp_build))
        non_updated_rs_that_should_exist_in_EVA = []
    

def ensure_new_rs_not_in_eva_human_accession_db(mongo_connection_handle: MongoClient, eva_production_human_dbsnp_build,
                                                    rs_id=None):
    global updated_rs_that_should_not_exist_in_EVA, min_rs_records_to_check

    num_rs_records_to_check = min_rs_records_to_check
    if rs_id:
        updated_rs_that_should_not_exist_in_EVA.append(rs_id)
    else:
        num_rs_records_to_check = len(updated_rs_that_should_not_exist_in_EVA)
    if len(updated_rs_that_should_not_exist_in_EVA) == num_rs_records_to_check or \
            (rs_id is None and num_rs_records_to_check > 0):
        results = get_mongo_query_result(updated_rs_that_should_not_exist_in_EVA, mongo_connection_handle)
        if len(results) > 0:
            logger.error("Found RS IDs {0} in EVA production even though they were marked in JSON "
                         "as released after dbSNP build {1} currently in EVA production!!"
                         .format(results, eva_production_human_dbsnp_build))
        updated_rs_that_should_not_exist_in_EVA = []


def is_rs_id_mapped_to_assembly(rs_record, eva_production_human_dbsnp_assembly):
    try:
        for placement in rs_record["primary_snapshot_data"]["placements_with_allele"]:
            if placement["is_ptlp"]:
                for assembly_info in placement["placement_annot"]["seq_id_traits_by_assembly"]:
                    if assembly_info["assembly_accession"] == eva_production_human_dbsnp_assembly:
                        return True
    except KeyError as ex:
        logger.error(traceback.format_exc())
        return False
    return False


def check_RS_release_JSON_assumptions(private_config_xml_file, release_json_file, eva_production_human_dbsnp_build,
                                      eva_production_human_dbsnp_assembly):
    with bz2.open(release_json_file) as release_json_file_handle, \
            MongoClient(get_mongo_uri_for_eva_profile("production", private_config_xml_file)) \
                    as mongo_connection_handle:
        line_index = 0
        for json_line in release_json_file_handle:
            if line_index % 100000 == 0:
                logger.info("Processed {0} records...".format(line_index))
            line_index += 1
            rs_record = json.loads(json_line.decode("utf-8").strip())
            if not (is_rs_id_mapped_to_assembly(rs_record, eva_production_human_dbsnp_assembly)):
                logger.error("RS ID {0} is not mapped to assembly {1}".format(rs_record["refsnp_id"],
                                                                              eva_production_human_dbsnp_assembly))
                continue
            rs_id = int(rs_record["refsnp_id"])

            if "support" in rs_record["primary_snapshot_data"]:
                support_record_count = 0
                support_record_last_updated_builds = set()
                for support_record in rs_record["present_obs_movements"]:
                    if "last_added_to_this_rs" in support_record:
                        support_record_count += 1
                        support_record_last_updated_builds.add(int(support_record["last_added_to_this_rs"]))
                if support_record_count == 0:
                    logger.error("Support record not found for RS {0} at line {1}!!".format(rs_id, line_index))
                    continue
                rs_last_updated_build = min(support_record_last_updated_builds)
                # Ensure newly added RS are not in production
                if rs_last_updated_build > eva_production_human_dbsnp_build:
                    ensure_new_rs_not_in_eva_human_accession_db(mongo_connection_handle,
                                                                    eva_production_human_dbsnp_build, rs_id)
                # Ensure RS in previous builds are present in production
                if rs_last_updated_build <= eva_production_human_dbsnp_build:
                    ensure_existing_rs_in_eva_human_accession_db(mongo_connection_handle, eva_production_human_dbsnp_build,
                                                                 rs_id)
        
        ensure_new_rs_not_in_eva_human_accession_db(mongo_connection_handle, eva_production_human_dbsnp_build)
        ensure_existing_rs_in_eva_human_accession_db(mongo_connection_handle, eva_production_human_dbsnp_build)


@click.option("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
@click.option("--release-json-file", help="ex: /path/to/release/json_file.json", required=True)
@click.option("--eva-production-human-dbsnp-build",
              help="Most recent dbSNP human release build in EVA production (ex: 152)", type=int, required=True)
@click.option("--eva-production-human-dbsnp-assembly",
              help="Most recent dbSNP human assembly in EVA production (ex: GCF_000001405.38)", required=True)
@click.command()
def main(private_config_xml_file, release_json_file, eva_production_human_dbsnp_build,
         eva_production_human_dbsnp_assembly):
    check_RS_release_JSON_assumptions(private_config_xml_file, release_json_file, eva_production_human_dbsnp_build,
                                      eva_production_human_dbsnp_assembly)


if __name__ == "__main__":
    main()
