import argparse
import copy
import hashlib
import os
import traceback
from datetime import datetime
from itertools import islice

import pymongo.read_preferences
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern
from pymongo.write_concern import WriteConcern

logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)

DBSNP_CLUSTERED_VARIANT_ENTITY = 'dbsnpClusteredVariantEntity'
EVA_CLUSTERED_VARIANT_ENTITY = 'clusteredVariantEntity'
DBSNP_SUBMITTED_VARIANT_ENTITY = 'dbsnpSubmittedVariantEntity'
EVA_SUBMITTED_VARIANT_ENTITY = 'submittedVariantEntity'
DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY = 'dbsnpClusteredVariantOperationEntity'
EVA_CLUSTERED_VARIANT_OPERATION_ENTITY = 'clusteredVariantOperationEntity'


def get_SHA1(text):
    h = hashlib.sha1()
    h.update(text.encode())
    return h.hexdigest().upper()


def get_submitted_SHA1(variant_rec):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    keys = ['seq', 'study', 'contig', 'start', 'ref', 'alt']
    return get_SHA1('_'.join([str(variant_rec[key]) for key in keys]))


def get_clustered_SHA1(variant_rec):
    """Calculate the SHA1 digest from the asm, contig, start and type attributes of the variant"""
    keys = ['asm', 'contig', 'start', 'type']
    return get_SHA1('_'.join([str(variant_rec[key]) for key in keys]))


def fix_discordant_variants(mongo_source, assembly, rs_file, batch_size=1000):
    logger.info(f"\n\nStarted processing assembly : {assembly}")

    with open(rs_file, 'r') as rs_file:
        while True:
            rs_list_to_process = [int(rs.strip()) for rs in list(islice(rs_file, batch_size))]
            if not rs_list_to_process:
                break

            rs_list = rs_list_to_process.copy()
            rs_list_to_process.clear()

            while rs_list:
                logger.info(f"Processing Batch. Num_of_RS in batch : {len(rs_list)} \nRS_list -> {rs_list}")

                all_rs_variants = get_rs_variants(mongo_source, assembly, rs_list)
                dbsnp_ss_variants, eva_ss_variants, all_ss_variants = get_ss_variants(mongo_source, assembly, rs_list)
                all_events = {}

                for rs in rs_list:
                    logger.info(f"Started Processing RS {rs}")

                    if rs not in all_rs_variants or not all_rs_variants[rs]:
                        logger.error(f"No RS variant could be found for RS {rs}")

                        if not all_events:
                            all_events = get_rs_events(mongo_source, rs_list, assembly)

                        # If no RS found in DB, check if the original RS has been merged to some other RS,
                        # if yes, add the new RS to the list for processing
                        if rs not in all_events:
                            logger.error(f"No RS merge events could be found for RS {rs}")
                            continue

                        rs_events = all_events[rs]
                        for event in rs_events:
                            if event['accession'] == rs and event['eventType'] == 'MERGED':
                                merge_into_RS = event['mergeInto']
                                logger.error(
                                    f"RS {rs} has been merged into RS {merge_into_RS}. Adding {merge_into_RS} to the list")
                                rs_list_to_process.append(merge_into_RS)
                                break
                        continue

                    rs_records = all_rs_variants[rs]

                    rs_without_map_weight = get_rs_without_map_weight(rs_records)
                    if not rs_without_map_weight:
                        logger.error(f"All variants for {rs} are map-weighted : \n{rs_records}")
                        continue
                    elif len(rs_without_map_weight) > 1:
                        logger.error(f"More than one variant without map-weight found for RS {rs} :"
                                     f"\n {rs_without_map_weight}")
                        continue

                    if rs not in all_ss_variants or not all_ss_variants[rs]:
                        logger.error(f"No original SS record found for RS {rs}")
                        continue

                    rs_variant = rs_without_map_weight[0]
                    ss_records = all_ss_variants[rs]

                    if ss_has_map_weight(ss_records):
                        logger.error(f"RS {rs} has SS with map weight.{ss_records}")
                        continue

                    if check_all_ss_has_same_info(ss_records):
                        if rs_variant['start'] == ss_records[0]['start']:
                            logger.error(f"RS {rs} and original SS's Start matches. Nothing to do")
                            continue

                        logger.info(f"Correct Discordant variants for RS {rs}")
                        merged_rs, merged_into = correct_discordant_rs_and_insert_into_db(rs_variant, ss_records,
                                                                                         assembly)

                        # check if any merge has happened and
                        # if the merged_rs and the merged_into rs are both in the same batch
                        # if not we don't need to do anything
                        # if yes
                        if merged_rs is not None and merged_rs in rs_list and merged_into in rs_list:
                            # check if the merged_rs is present later in the list (after current rs),
                            # if yes remove it and put the RS it merged into in next batch for processing
                            if rs_list.index(merged_rs) > rs_list.index(rs):
                                logger.info(f"Removing rs {merged_rs} from current batch as it is merged "
                                            f"and putting rs {merged_into} for processing in the next batch")
                                rs_list.remove(merged_rs)
                                rs_list_to_process.append(merged_into)
                            # if the current RS is the one being merged and the RS it got merged into is in the same
                            # batch later, we need to update the info associated with that merged_into RS in memory
                            # as some of it might have changed because of the merge
                            # The easiest to do this is just put it into the next batch
                            else:
                                logger.info(f"Removing rs {merged_into} from current batch "
                                            f"as it is involved in merge event and putting it in the next batch")
                                rs_list.remove(merged_into)
                                rs_list_to_process.append(merged_into)

                    else:
                        logger.error(
                            f"For RS {rs}, Not all original SS has same info. Case for Split: \nSS Records {ss_records}")

                rs_list = rs_list_to_process.copy()
                rs_list_to_process.clear()


def correct_discordant_rs_and_insert_into_db(rs_variant, ss_records, assembly):
    rs_with_new_start = copy.copy(rs_variant)
    rs_with_new_start['start'] = ss_records[0]['start']
    rs_with_new_start['_id'] = get_clustered_SHA1(rs_with_new_start)

    variant_in_dbsnp, variant_in_eva = check_for_hash_collision(rs_with_new_start['_id'])

    if variant_in_dbsnp or variant_in_eva:
        logger.warn(f"Hash collision will occur for RS {rs_variant['accession']} "
                    f"with RS {variant_in_dbsnp['accession'] if variant_in_dbsnp else variant_in_eva['accession']}")
        return resolve_collision_and_insert_rs(rs_variant, rs_with_new_start, variant_in_dbsnp, variant_in_eva,
                                               assembly)
    else:
        logger.info(f"No hash collision for RS {rs_variant['accession']}")
        dbsnp_cve_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_CLUSTERED_VARIANT_ENTITY]
        # delete original rs
        logger.info(f"delete rs with wrong start : {rs_variant}")
        dbsnp_cve_collection.with_options(write_concern=WriteConcern("majority")) \
            .delete_one({'_id': rs_variant['_id']})
        # insert rs with new
        logger.info(f"Insert rs with new start and id(hash): {rs_with_new_start}")
        dbsnp_cve_collection.with_options(write_concern=WriteConcern("majority")) \
            .insert_one(rs_with_new_start)

        return None, None


def resolve_collision_and_insert_rs(rs_variant, rs_with_new_start, variant_in_dbsnp, variant_in_eva, assembly):
    dbsnp_cve_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_CLUSTERED_VARIANT_ENTITY]
    dbsnp_cvoe_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY]
    eva_cve_collection = mongo_source.mongo_handle[mongo_source.db_name][EVA_CLUSTERED_VARIANT_ENTITY]
    eva_cvoe_collection = mongo_source.mongo_handle[mongo_source.db_name][EVA_CLUSTERED_VARIANT_OPERATION_ENTITY]

    variant_in_db = variant_in_dbsnp if variant_in_dbsnp else variant_in_eva

    # For priority refer to:
    # https://github.com/EBIvariation/eva-accession/blob/0b2ae4cdb6f74152c5443c3831c02c1d76cf93f9/eva-accession-clustering/src/main/java/uk/ac/ebi/eva/accession/clustering/batch/io/ClusteredVariantMergingPolicy.java#L40
    if rs_with_new_start['accession'] < variant_in_db['accession']:
        if variant_in_dbsnp:
            merge_event = create_merge_event(variant_in_db, rs_with_new_start)
            dbsnp_cvoe_collection.with_options(write_concern=WriteConcern("majority")).insert_one(merge_event)
            logger.info(
                f"delete rs with wrong start and the one being merged: \nWrong start: {rs_variant} \nMerged: {variant_in_db}")
            dbsnp_cve_collection.with_options(write_concern=WriteConcern("majority")) \
                .delete_many({'_id': {'$in': [rs_variant['_id'], variant_in_db['_id']]}})
        else:
            merge_event = create_merge_event(variant_in_db, rs_with_new_start)
            eva_cvoe_collection.with_options(write_concern=WriteConcern("majority")).insert_one(merge_event)
            logger.info(
                f"delete rs with wrong start and the one being merged: \nWrong start: {rs_variant} \nMerged: {variant_in_db}")
            dbsnp_cve_collection.with_options(write_concern=WriteConcern("majority")) \
                .delete_many({'_id': {'$in': [rs_variant['_id']]}})
            eva_cve_collection.with_options(write_concern=WriteConcern("majority")) \
                .delete_many({'_id': {'$in': [variant_in_db['_id']]}})

        logger.info(f"insert rs with new start and hash : {rs_with_new_start}")
        dbsnp_cve_collection.with_options(write_concern=WriteConcern("majority")).insert_one(rs_with_new_start)

        update_ss_with_new_rs(variant_in_db['accession'], rs_with_new_start['accession'], assembly)

        return variant_in_db['accession'], rs_with_new_start['accession']

    else:
        merge_event = create_merge_event(rs_with_new_start, variant_in_db)
        dbsnp_cvoe_collection.with_options(write_concern=WriteConcern("majority")).insert_one(merge_event)

        logger.info(f"delete rs with wrong start: {rs_variant}")
        dbsnp_cve_collection.with_options(write_concern=WriteConcern("majority")).delete_one({'_id': rs_variant['_id']})

        update_ss_with_new_rs(rs_with_new_start['accession'], variant_in_db['accession'], assembly)

        return rs_with_new_start['accession'], variant_in_db['accession']


def create_merge_event(variant_merged, variant_retained):
    logger.info(
        f"creating merge event for accession: {variant_merged['accession']} mergeInto: {variant_retained['accession']}")
    merge_event = {
        "_id": f"EVA2850_MERGED_{variant_merged['accession']}_{variant_retained['_id']}",
        "eventType": "MERGED",
        "accession": variant_merged['accession'],
        "mergeInto": variant_retained['accession'],
        "reason": "EVA2850: Merged because of RS discordant correction",
        "inactiveObjects": [variant_merged],
        "createdDate": datetime.now()
    }

    return merge_event


def update_ss_with_new_rs(old_rs, new_rs, assembly):
    dbsnp_sve_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_SUBMITTED_VARIANT_ENTITY]
    eva_sve_collection = mongo_source.mongo_handle[mongo_source.db_name][EVA_SUBMITTED_VARIANT_ENTITY]

    logger.info(f"updating submittedVariantEntity with old_rs: {old_rs} to new_rs: {new_rs}")

    filter_query = {'rs': old_rs, 'seq': assembly}
    update_value = {'$set': {'rs': new_rs}}

    dbsnp_sve_collection.with_options(read_concern=ReadConcern("majority"),
                                      read_preference=pymongo.ReadPreference.PRIMARY,
                                      write_concern=WriteConcern("majority")) \
        .update_many(filter_query, update_value)
    eva_sve_collection.with_options(read_concern=ReadConcern("majority"),
                                    read_preference=pymongo.ReadPreference.PRIMARY,
                                    write_concern=WriteConcern("majority")) \
        .update_many(filter_query, update_value)


def get_rs_variants(mongo_source, assembly, rs_list):
    rs_filter_criteria = {'asm': assembly, 'accession': {'$in': rs_list}}
    rs_variants = get_variants(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY, rs_filter_criteria, 'accession')
    return rs_variants


def get_ss_variants(mongo_source, assembly, rs_list):
    ss_filter_criteria = {'seq': assembly, 'rs': {'$in': rs_list}, '$or': [{"allelesMatch": {"$exists": False}},
                                                                           {"allelesMatch": True}]}
    dbsnp_ss_variants = get_variants(mongo_source, DBSNP_SUBMITTED_VARIANT_ENTITY, ss_filter_criteria, 'rs')
    eva_ss_variants = get_variants(mongo_source, EVA_SUBMITTED_VARIANT_ENTITY, ss_filter_criteria, 'rs')
    all_ss_variants = merge_all_records(dbsnp_ss_variants, eva_ss_variants)
    return dbsnp_ss_variants, eva_ss_variants, all_ss_variants


def get_rs_events(mongo_source, rs_list, asm):
    event_filter_criteria = {"inactiveObjects.asm": asm, 'accession': {'$in': rs_list}, 'eventType': 'MERGED'}
    rs_events = get_events(mongo_source, DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY, event_filter_criteria)
    return rs_events


def get_events(mongo_source, collection_name, filter_criteria):
    records = {}
    try:
        for event in find_documents(mongo_source, collection_name, filter_criteria):
            if event['accession'] in records:
                records[event['accession']].append(event)
            else:
                records[event['accession']] = [event]
            if 'mergeInto' in event and event['mergeInto'] in records:
                records[event['mergeInto']].append(event)
            elif 'mergeInto' in event and event['mergeInto'] not in records:
                records[event['mergeInto']] = [event]
            if 'splitInto' in event and event['splitInto'] in records:
                records[event['splitInto']].append(event)
            elif 'splitInto' in event and event['splitInto'] not in records:
                records[event['splitInto']] = [event]
    except Exception as e:
        print(traceback.format_exc())
        raise e

    return records


def get_variants(mongo_source, collection_name, filter_criteria, key):
    records = {}
    try:
        for variant in find_documents(mongo_source, collection_name, filter_criteria):
            if variant[key] in records:
                records[variant[key]].append(variant)
            else:
                records[variant[key]] = [variant]
    except Exception as e:
        print(traceback.format_exc())
        raise e

    return records


def find_documents(mongo_source, collection_name, filter_criteria):
    collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
    cursor = collection.with_options(read_concern=ReadConcern("majority"),
                                     read_preference=pymongo.ReadPreference.PRIMARY) \
        .find(filter_criteria, no_cursor_timeout=True)
    records = []
    try:
        for result in cursor:
            records.append(result)
    except Exception as e:
        logger.exception(traceback.format_exc())
        raise e
    finally:
        cursor.close()

    return records


def merge_all_records(dbsnp_records, eva_records):
    all_keys = dbsnp_records.keys() | eva_records.keys()
    all_rs_records = {}
    for key in all_keys:
        all_rs_records[key] = []
        if key in dbsnp_records:
            all_rs_records[key].extend(dbsnp_records[key])
        if key in eva_records:
            all_rs_records[key].extend(eva_records[key])

    return all_rs_records


def get_rs_without_map_weight(variants_list):
    variants_without_map_weight = []
    for variant in variants_list:
        if 'mapWeight' not in variant:
            variants_without_map_weight.append(variant)
        elif variant['mapWeight'] is None:
            variants_without_map_weight.append(variant)

    return variants_without_map_weight


def ss_has_map_weight(ss_list):
    for ss in ss_list:
        if 'mapWeight' in ss:
            return True

    return False


def check_all_ss_has_same_info(ss_list):
    contig = set([ss['contig'] for ss in ss_list])
    start = set([ss['start'] for ss in ss_list])

    if len(contig) > 1 or len(start) > 1:
        return False
    else:
        return True


def check_for_hash_collision(id):
    hash_collision_filter_criteria = {'_id': id}
    dbsnp_collision_records = find_documents(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY,
                                             hash_collision_filter_criteria)
    eva_collision_records = find_documents(mongo_source, EVA_CLUSTERED_VARIANT_ENTITY, hash_collision_filter_criteria)

    if dbsnp_collision_records and not eva_collision_records:
        return dbsnp_collision_records[0], None
    elif not dbsnp_collision_records and eva_collision_records:
        return None, eva_collision_records[0]
    elif not dbsnp_collision_records and not eva_collision_records:
        return None, None
    else:
        raise Exception(f"CVE variant with {id}  is present in both dbsnp cve and eva cve")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find discordant variants', add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)",
                        required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--discordant-rs-dir", help="File containing discordant rs ids", required=True)
    args = parser.parse_args()

    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva2959_accession_sharded")

    # Sometimes the cluster creates hidden files like .nfs<numeric ID> in the folder when a file is being read. So we just want to focus on the ones that start with GCA
    all_files = [os.path.join(args.discordant_rs_dir, filename) for filename in os.listdir(args.discordant_rs_dir) if filename.startswith("GCA")]
    for file in sorted(all_files, key=lambda x: os.stat(x).st_size):
        assembly = os.path.basename(file)
        fix_discordant_variants(mongo_source, assembly, os.path.join(args.discordant_rs_dir, assembly))

    logger.info(f"Process Finished")
