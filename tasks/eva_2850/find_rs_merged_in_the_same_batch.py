import argparse
import os
import signal
from collections import defaultdict

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from ebi_eva_common_pyutils.network_utils import forward_remote_port_to_local_port, get_available_local_port
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query

from tasks.eva_2850.fix_discordant_variants import get_variants, DBSNP_SUBMITTED_VARIANT_ENTITY, \
    EVA_SUBMITTED_VARIANT_ENTITY, merge_all_records, DBSNP_CLUSTERED_VARIANT_ENTITY, get_SHA1, \
    DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY, find_documents

logger = logging_config.get_logger(__name__)


def merged_rs_ids_present_in_same_batch(all_log_files):
    merged_rs_ids = {}
    merged_rs_with_further_hash_collision = []

    for log_file_path in sorted(all_log_files, key=lambda x: os.stat(x).st_size):
        with open(log_file_path, 'r') as log_file:
            curr_rs_id_batch = []
            for line in log_file:
                if "RS_list ->" in line:
                    rs_id_list = line[line.find('['):].split(",")
                    rs_id_list = [s.replace("[", "") for s in rs_id_list]
                    rs_id_list = [s.replace("]", "") for s in rs_id_list]
                    rs_id_list = [s.strip() for s in rs_id_list]
                    # only one rs in batch
                    if len(rs_id_list) == 1:
                        continue
                    curr_rs_id_batch = rs_id_list
                if "creating merge event for" in line:
                    merged_rs = line[line.index("creating merge event for"):].split(":")[1].split(" ")[1].strip()
                    merged_into = line[line.index("creating merge event for"):].split(":")[2].strip()
                    # check if the variants are involved in another hash collision
                    if merged_rs in merged_rs_ids:
                        merged_rs_with_further_hash_collision.append(merged_rs)
                    if merged_into in merged_rs_ids:
                        merged_rs_with_further_hash_collision.append(merged_into)

                    if merged_rs in curr_rs_id_batch and merged_into in curr_rs_id_batch:
                        merged_rs_ids[merged_rs] = merged_into

    logger.info(f"All rs ids merged in same batch. Total=> {len(merged_rs_ids.keys())} RS=> {merged_rs_ids}")
    rs_not_involved_in_further_merge = list(set(merged_rs_ids.keys()) - set(merged_rs_with_further_hash_collision))
    logger.info(f"RS ids not involved in further merge events: {rs_not_involved_in_further_merge}")
    logger.info(f"RS ids involved in further merge events: {set(merged_rs_with_further_hash_collision)}")

    # get logs for each rs involved
    rs_logs = defaultdict(list)
    for log_file_path in sorted(all_log_files, key=lambda x: os.stat(x).st_size):
        with open(log_file_path, 'r') as log_file:
            for line in log_file:
                split_line = [s.strip().replace(",", "") for s in line.split(" ")]
                if "RS_list ->" in line:
                    continue
                for rs in merged_rs_ids:
                    if rs in split_line or merged_rs_ids[rs] in split_line:
                        rs_logs[rs].append(line)

    # print logs for each rs
    for rs, logs in rs_logs.items():
        print(f"------------------------------------------------------------------------------------------------------")
        print(f"logs grep command : grep -n -E \"{rs}|{merged_rs_ids[rs]}\"\n")
        for log_line in logs:
            print(log_line.replace("\n", ""))


def correct_sve_with_wrong_rs(all_log_files, mongo_source, private_config_xml_file):
    rs_list = []

    # go through each of the log files and find all rs for which sve has been updated
    for log_file_path in sorted(all_log_files, key=lambda x: os.stat(x).st_size):
        with open(log_file_path, 'r') as log_file:
            for line in log_file:
                # get the new_rs from the line where sve with old_rs updated to new_rs
                if "updating submittedVariantEntity with old_rs:" in line:
                    # old_rs = line[line.index("updating submittedVariantEntity with old_rs:"):].split(":")[1].strip() \
                    #     .split(" ")[0].strip()
                    new_rs = int(
                        line[line.index("updating submittedVariantEntity with old_rs:"):].split(":")[2].strip())
                    rs_list.append(new_rs)

    all_rs_variants = get_rs_variants(mongo_source, list(set(rs_list)))
    rs_not_found_in_db = []
    for rs in rs_list:
        # rs might have been merged  or deleted because of some other remediation
        if rs not in all_rs_variants:
            rs_not_found_in_db.append(rs)
    logger.info(f"RS not found in db: Total => {len(rs_not_found_in_db)} RS=> {rs_not_found_in_db}")
    logger.info(f"Checking if these have been merged to some other RS, if yes add those RS to the list")
    merge_events = get_rs_merge_events(mongo_source, rs_not_found_in_db)

    rs_not_found_in_db_not_merged = []
    for rs in rs_not_found_in_db:
        if rs not in merge_events:
            rs_not_found_in_db_not_merged.append(rs)
        else:
            rs_list.append(merge_events[rs][0]['mergeInto'])

    logger.info(f"RS which was not found in db and could not find in any merge event: "
                f"Total => {len(rs_not_found_in_db_not_merged)}, RS List => {rs_not_found_in_db_not_merged}")

    _, _, all_ss_variants = get_ss_variants(mongo_source, list(set(rs_list)))
    rs_for_which_no_sve_found = []
    for rs in rs_list:
        # no sve found for rs
        if rs not in all_ss_variants:
            rs_for_which_no_sve_found.append(rs)

    logger.info(f"No SVE found for following RS: {rs_for_which_no_sve_found}")

    print("\n")
    ss_list = [sve for sve_list in all_ss_variants.values() for sve in sve_list]
    logger.info(f"No of ss needed to be checked for correction : {len(ss_list)}")

    all_hashes_from_sve = list(set([get_rs_hash_from_sve(sve) for sve in ss_list]))
    all_cve_for_sve = get_rs_variants_with_hashes(mongo_source, all_hashes_from_sve, DBSNP_CLUSTERED_VARIANT_ENTITY)

    sve_for_which_cve_not_found = check_sve_with_cve_for_correction(ss_list, all_cve_for_sve)

    if sve_for_which_cve_not_found:
        logger.info("Getting SVE from tempmongo")
        all_sve_from_tempmongo = get_sve_from_tempmongo(sve_for_which_cve_not_found, private_config_xml_file)
        logger.info(f"Total SVE found in tempmongo: {len(all_sve_from_tempmongo)}")
        sve_not_found_in_tempmongo = check_sve_with_sve_for_correction(sve_for_which_cve_not_found,
                                                                       all_sve_from_tempmongo)

        if sve_not_found_in_tempmongo:
            logger.info(f"There are some sve which could not be found in tempmongo")
            logger.info(f"total => {len(sve_not_found_in_tempmongo)}")
            for sve in sve_not_found_in_tempmongo:
                logger.info(f"{sve}")

    logger.info("Finished")


def get_sve_from_tempmongo(sve_list_for_tempmongo, private_config_xml_file):
    tax_sve_list = defaultdict(list)
    for sve in sve_list_for_tempmongo:
        tax_sve_list[sve['tax']].append(sve)

    tax_tempmongo = get_tax_tempmongo(private_config_xml_file)

    sve_from_tempmongo = {}
    for tax, sve_list_for_tax in tax_sve_list.items():
        for tempmongo in tax_tempmongo[tax]:
            logger.info(
                    f"Connecting to tempmongo {tempmongo} for taxonomy {tax} with sve list of len {len(sve_list_for_tax)}")
            all_sve_ids_for_tempmongo = list(set([sve['_id'] for sve in sve_list_for_tax]))
            sve_from_tempmongo.update(get_sve_for_taxonomy_from_tempmongo(all_sve_ids_for_tempmongo, tax, tempmongo))


    return sve_from_tempmongo


def get_sve_for_taxonomy_from_tempmongo(all_sve_ids_for_tempmongo, taxonomy, tempmongo_instance):
    MONGO_PORT = 27017
    local_forwarded_port = get_available_local_port(MONGO_PORT)
    try:
        logger.info("Forwarding remote MongoDB port 27017 to local port {0}...".format(local_forwarded_port))
        port_forwarding_process_id = forward_remote_port_to_local_port(tempmongo_instance, MONGO_PORT,
                                                                       local_forwarded_port)
        mongo_source = MongoDatabase(uri=f"mongodb://localhost:{local_forwarded_port}/?authSource=admin",
                                     secrets_file=None, db_name=f"acc_{taxonomy}")
        dbsnp_sve_from_tempmongo = get_sve_variants_with_hashes(mongo_source, all_sve_ids_for_tempmongo,
                                                                DBSNP_SUBMITTED_VARIANT_ENTITY)
        eva_sve_from_tempmongo = get_sve_variants_with_hashes(mongo_source, all_sve_ids_for_tempmongo,
                                                              EVA_SUBMITTED_VARIANT_ENTITY)
        all_sve_from_tempmongo = merge_all_records(dbsnp_sve_from_tempmongo, eva_sve_from_tempmongo)
        logger.info(
            f"For {taxonomy} in {tempmongo_instance} retrieved a total of  {len(all_sve_from_tempmongo)} SVE")
        return all_sve_from_tempmongo
    finally:
        close_mongo_port_to_tempmongo(port_forwarding_process_id)


def close_mongo_port_to_tempmongo(port_forwarding_process_id):
    os.kill(port_forwarding_process_id, signal.SIGTERM)
    os.system('echo -e "Killed port forwarding from remote port with signal 1 - SIGTERM. '
              '\\033[31;1;4mIGNORE OS MESSAGE '  # escape sequences for bold red and underlined text
              '\'Killed by Signal 1\' in the preceding/following text\\033[0m".')


def get_tax_tempmongo(private_config_xml_file):
    tax_tempmongo = defaultdict(list)
    with get_metadata_connection_handle('production_processing', private_config_xml_file) as pg_conn:
        query = f"select distinct taxonomy, tempmongo_instance from eva_progress_tracker.clustering_release_tracker crt " \
                f"where release_version = 3 and tempmongo_instance is not null order by taxonomy  "
        for taxonomy, tempmongo in get_all_results_for_query(pg_conn, query):
            tax_tempmongo[taxonomy].append(tempmongo)

    return tax_tempmongo


def check_sve_with_cve_for_correction(ss_list, all_cve_for_sve):
    sve_for_which_cve_not_found = []
    sve_with_wrong_rs = []

    for sve in ss_list:
        cve_hash_from_sve = get_rs_hash_from_sve(sve)
        if cve_hash_from_sve not in all_cve_for_sve:
            sve_for_which_cve_not_found.append(sve)
        else:
            cve = all_cve_for_sve[cve_hash_from_sve]
            if sve['rs'] != cve[0]['accession']:
                sve_with_wrong_rs.append(sve)
                logger.info(f"SVE: {sve} \nCVE: {cve[0]}")
                logger.info(f"SVE RS {sve['rs']} does not match with CVE accession {cve[0]['accession']}")

    logger.info(f"SVE needs to be corrected: {len(sve_with_wrong_rs)}")
    logger.info(f"SVE for which CVE is not found: {len(sve_for_which_cve_not_found)}")

    return sve_for_which_cve_not_found


def check_sve_with_sve_for_correction(ss_list, all_sve_from_tempmongo):
    sve_not_found = []
    sve_with_wrong_rs = []

    for sve in ss_list:
        if sve['_id'] not in all_sve_from_tempmongo:
            sve_not_found.append(sve)
        else:
            sve_from_tempmongo = all_sve_from_tempmongo[sve['_id']]
            if sve['rs'] != sve_from_tempmongo[0]['rs']:
                sve_with_wrong_rs.append(sve)
                logger.info(f"SVE: {sve} \nSVE from Tempmongo: {sve_from_tempmongo[0]}")
                logger.info(
                    f"SVE RS {sve['rs']} does not match with SVE RS from tempmongo {sve_from_tempmongo[0]['rs']}")

    logger.info(f"SVE needs to be corrected: {len(sve_with_wrong_rs)}")
    logger.info(f"SVE for which CVE is not found: {len(sve_not_found)}")

    return sve_not_found


def check_if_any_newly_added_rs_has_collision_in_eva_cve(all_log_files, mongo_source):
    inserted_rs_id_list = []
    for log_file_path in sorted(all_log_files, key=lambda x: os.stat(x).st_size):
        with open(log_file_path, 'r') as log_file:
            for line in log_file:
                if "Insert rs with new start and id" in line:
                    id_field = line[line.index("{"):].split(",")[0]
                    id = id_field.split(":")[1].replace("'", "").strip()
                    inserted_rs_id_list.append(id)
                    continue
                if "insert rs with new start and hash :" in line:
                    id_field = line[line.index("{"):].split(",")[0]
                    id = id_field.split(":")[1].replace("'", "").strip()
                    inserted_rs_id_list.append(id)

    eva_rs_variants = get_rs_variants_with_hashes(mongo_source, inserted_rs_id_list, "clusteredVariantEntity")
    logger.info(f"No of rs ids found in eva : {len(eva_rs_variants)}")
    for id, cve in eva_rs_variants.items():
        logger.info(f"{cve}")


def check_if_all_processed_rs_was_supposed_to_be_processed(all_log_files, mongo_source):
    asm_rs_list = {}
    curr_asm = ""
    for log_file_path in sorted(all_log_files, key=lambda x: os.stat(x).st_size):
        with open(log_file_path, 'r') as log_file:
            for line in log_file:
                if 'Started processing assembly :' in line:
                    curr_asm = line.split(':')[1].strip()
                    asm_rs_list[curr_asm] = []
                if "Correct Discordant variants for RS" in line:
                    rs = line[line.index("Correct Discordant variants for RS"):].split(" ")[-1].strip()
                    asm_rs_list[curr_asm].append(int(rs))

    for asm, rs_list in asm_rs_list.items():
        rs_variants = get_rs_variants_with_asm(mongo_source, asm, list(set(rs_list)))
        for rs, cve_list in rs_variants.items():
            if len(cve_list) > 1:
                logger.info(f"RS {rs}, cve_list_length : {len(cve_list)}")
                for cve in cve_list:
                    logger.info(f"{cve}")

    logger.info("finished")


def get_rs_with_hash(rs_list, hash):
    for cve in rs_list:
        if cve['_id'] == hash:
            return cve

    return None


def get_rs_hash_from_sve(sve):
    asm = sve['seq']
    contig = sve['contig']
    start = sve['start']
    type = get_variant_type(sve)

    fields = [asm, contig, start, type]
    return get_SHA1('_'.join(str(field) for field in fields))


def get_variant_type(sve):
    if sve['ref'] and sve['alt']:
        if len(sve['ref']) == len(sve['alt']):
            if len(sve['ref']) == 1:
                return "SNV"
            else:
                return "MNV"
        else:
            return "INDEL"

    if sve['ref'] and not sve['alt']:
        return "DEL"

    if sve['alt'] and not sve['ref']:
        return "INS"


def get_ss_variants(mongo_source, rs_list):
    ss_filter_criteria = {'rs': {'$in': rs_list}}
    dbsnp_ss_variants = get_variants(mongo_source, DBSNP_SUBMITTED_VARIANT_ENTITY, ss_filter_criteria, 'rs')
    eva_ss_variants = get_variants(mongo_source, EVA_SUBMITTED_VARIANT_ENTITY, ss_filter_criteria, 'rs')
    all_ss_variants = merge_all_records(dbsnp_ss_variants, eva_ss_variants)
    return dbsnp_ss_variants, eva_ss_variants, all_ss_variants


def get_rs_variants_with_asm(mongo_source, assembly, rs_list):
    rs_filter_criteria = {'asm': assembly, 'accession': {'$in': rs_list}}
    rs_variants = get_variants(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY, rs_filter_criteria, 'accession')
    return rs_variants


def get_rs_variants(mongo_source, rs_list):
    rs_filter_criteria = {'accession': {'$in': rs_list}}
    rs_variants = get_variants(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY, rs_filter_criteria, 'accession')
    return rs_variants


def get_rs_merge_events(mongo_source, rs_list):
    event_filter_criteria = {'accession': {'$in': rs_list}, 'eventType': 'MERGED'}
    rs_events = get_events(mongo_source, DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY, event_filter_criteria)
    return rs_events


def get_events(mongo_source, collection_name, filter_criteria):
    records = {}
    for event in find_documents(mongo_source, collection_name, filter_criteria):
        if event['accession'] in records:
            records[event['accession']].append(event)
        else:
            records[event['accession']] = [event]

    return records


def get_rs_variants_with_hashes(mongo_source, rs_hash_list, collection):
    rs_filter_criteria = {'_id': {'$in': rs_hash_list}}
    rs_variants = get_variants(mongo_source, collection, rs_filter_criteria, '_id')
    return rs_variants


def get_sve_variants_with_hashes(mongo_source, sve_id_list, collection):
    sve_filter_criteria = {'_id': {'$in': sve_id_list}}
    sve_variants = get_variants(mongo_source, collection, sve_filter_criteria, '_id')
    return sve_variants


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Find variants which are part of the same batch and involved in merge', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--log-file-dir", help="File containing discordant rs ids", required=True)
    args = parser.parse_args()

    # there are 2 different log files
    all_log_files = [os.path.join(args.log_file_dir, filename) for filename in os.listdir(args.log_file_dir)]

    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")

    merged_rs_ids_present_in_same_batch(all_log_files)
    correct_sve_with_wrong_rs(all_log_files, mongo_source, args.private_config_xml_file)
    check_if_any_newly_added_rs_has_collision_in_eva_cve(all_log_files, mongo_source)
    check_if_all_processed_rs_was_supposed_to_be_processed(all_log_files, mongo_source)
