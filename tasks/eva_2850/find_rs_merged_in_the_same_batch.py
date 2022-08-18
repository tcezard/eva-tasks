import argparse
import os

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase

from tasks.eva_2850.fix_discordant_variants import get_variants, DBSNP_SUBMITTED_VARIANT_ENTITY, \
    EVA_SUBMITTED_VARIANT_ENTITY, merge_all_records, DBSNP_CLUSTERED_VARIANT_ENTITY, get_SHA1, \
    DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY, get_events

logger = logging_config.get_logger(__name__)


def merged_rs_ids_present_in_same_batch(all_log_files):
    merged_rs_ids = []
    merged_rs_ids_no_further_hash_collision = []
    merged_rs_with_further_hash_collision = []

    for log_file_path in sorted(all_log_files, key=lambda x: os.stat(x).st_size):
        with open(log_file_path, 'r') as log_file:
            curr_rs_id_batch = []
            for line in log_file:
                if "RS_list ->" in line:
                    # logger.info(curr_rs_id_batch)
                    rs_id_list = line[line.find('['):].split(",")
                    rs_id_list = [s.replace("[", "") for s in rs_id_list]
                    rs_id_list = [s.replace("]", "") for s in rs_id_list]
                    rs_id_list = [s.strip() for s in rs_id_list]
                    # only one rs in batch
                    if len(rs_id_list) == 1:
                        continue
                    curr_rs_id_batch = rs_id_list
                if "creating merge event for" in line:
                    # logger.info(line)
                    merged_rs = line[line.index("creating merge event for"):].split(":")[1].split(" ")[1].strip()
                    merged_into = line[line.index("creating merge event for"):].split(":")[2].strip()

                    # check if the variants (which should have been deleted) are involved in another hash collision
                    if merged_rs in merged_rs_ids:
                        merged_rs_with_further_hash_collision.append(merged_rs)
                    if merged_into in merged_rs_ids:
                        merged_rs_with_further_hash_collision.append(merged_into)

                    if merged_rs in curr_rs_id_batch and merged_into in curr_rs_id_batch:
                        # logger.info(curr_rs_id_batch)
                        if curr_rs_id_batch.index(merged_rs) > curr_rs_id_batch.index(merged_into):
                            logger.error(f"Merged RS : {merged_rs} Merged Into: {merged_into}")
                            # logger.error(f"grep -n -E \"{merged_rs}|{merged_into}\" {os.path.basename(log_file_path)}")
                            merged_rs_ids.append(merged_rs)

                if "No hash collision for RS" in line:
                    rs_id = line.split(" ")[-1].strip()
                    if rs_id in merged_rs_ids:
                        merged_rs_ids_no_further_hash_collision.append(rs_id)

    logger.info(f"All rs ids merged in same batch: {sorted(merged_rs_ids)}")
    logger.info(f"RS ids not involved in further merge events: "
                f"{sorted(set(merged_rs_ids_no_further_hash_collision) - set(merged_rs_with_further_hash_collision))}")
    logger.info(f"RS ids involved in further merge events: "
                f"{sorted(set(merged_rs_with_further_hash_collision) - set(merged_rs_ids_no_further_hash_collision))}")
    logger.info(f"RS ids processed twice - "
                f"once not involved in any merge and again involved in another merge event"
                f"{sorted(set(merged_rs_ids_no_further_hash_collision).intersection(set(merged_rs_with_further_hash_collision)))}")


def correct_sve_with_wrong_rs(all_log_files, mongo_source):
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
    sve_for_which_cve_not_found = check_ss_for_correction(ss_list, False)

    print("\n")
    logger.info(f"After start correction : (adding 1 to the start) ")
    print("\n")

    sve_for_which_cve_not_found = check_ss_for_correction(sve_for_which_cve_not_found, True)

    print("\nSVE for which CVE is not found even after 'start' correction: ")
    for sve in sve_for_which_cve_not_found:
        logger.info(sve)


def check_ss_for_correction(ss_list, start_correction):
    logger.info(f"No of ss in list : {len(ss_list)}")
    all_hashes_from_sve = list(set([get_rs_hash_from_sve(sve, start_correction) for sve in ss_list]))
    logger.info(f"No of unique cve hashes from sve:  {len(all_hashes_from_sve)}")
    all_cve_for_sve = get_rs_variants_with_hashes(mongo_source, all_hashes_from_sve)
    logger.info(f"No of cve found for sve using hashes:  {len(all_cve_for_sve)}")

    sve_for_which_cve_not_found = []
    sve_with_wrong_rs = []
    for sve in ss_list:
        cve_hash_from_sve = get_rs_hash_from_sve(sve, start_correction)
        if cve_hash_from_sve not in all_cve_for_sve:
            sve_for_which_cve_not_found.append(sve)
        else:
            cve = all_cve_for_sve[cve_hash_from_sve]
            if sve['rs'] != cve[0]['accession']:
                sve_with_wrong_rs.append(sve)
                # logger.info(f"SVE RS {sve['rs']} does not match with CVE accession {cve[0]['accession']}")

    logger.info(f"SVE needs to be corrected: {len(sve_with_wrong_rs)}")
    logger.info(f"SVE for which CVE is not found: {len(sve_for_which_cve_not_found)}")

    return sve_for_which_cve_not_found


def get_rs_with_hash(rs_list, hash):
    for cve in rs_list:
        if cve['_id'] == hash:
            return cve

    return None


def get_rs_hash_from_sve(sve, start_correction):
    asm = sve['seq']
    contig = sve['contig']
    if start_correction:
        start = sve['start'] + 1
    else:
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


def get_rs_variants(mongo_source, rs_list):
    rs_filter_criteria = {'accession': {'$in': rs_list}}
    rs_variants = get_variants(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY, rs_filter_criteria, 'accession')
    return rs_variants


def get_rs_merge_events(mongo_source, rs_list):
    event_filter_criteria = {'accession': {'$in': rs_list}}
    rs_events = get_events(mongo_source, DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY, event_filter_criteria)
    return rs_events


def get_rs_variants_with_hashes(mongo_source, rs_hash_list):
    rs_filter_criteria = {'_id': {'$in': rs_hash_list}}
    rs_variants = get_variants(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY, rs_filter_criteria, '_id')
    return rs_variants


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Find variants which are part of the same batch and involved in merge', add_help=False)
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

    # merged_rs_ids_present_in_same_batch(all_log_files)
    correct_sve_with_wrong_rs(all_log_files, mongo_source)
