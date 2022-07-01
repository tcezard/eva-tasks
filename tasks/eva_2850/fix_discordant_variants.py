import argparse
import copy
import hashlib
import os
import traceback
from itertools import islice

import pymongo
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern

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
            rs_list = [int(rs.strip()) for rs in list(islice(rs_file, batch_size))]
            if not rs_list:
                break

            logger.info(f"Processing Batch. Num_of_RS in batch : {len(rs_list)} \nRS_list -> {rs_list}")

            dbsnp_rs_variants, eva_rs_variants, all_rs_variants = get_rs_variants(mongo_source, assembly, rs_list)
            dbsnp_ss_variants, eva_ss_variants, all_ss_variants = get_ss_variants(mongo_source, assembly, rs_list)
            dbsnp_events, eva_events, all_events = get_rs_events(mongo_source, rs_list)

            for rs in rs_list:
                logger.info(f"Started Processing RS {rs}")

                if rs not in all_rs_variants:
                    logger.error(f"No RS variant could be found for RS {rs}")
                    rs_events = all_events[rs]
                    for event in rs_events:
                        if event['accession'] == rs and event['eventType'] == 'MERGED':
                            merge_into_RS = event['mergeInto']
                            logger.error(
                                f"RS {rs} has been merged into RS {merge_into_RS}. Adding {merge_into_RS} to the list")
                            rs_list.append(rs)
                            break
                    continue
                else:
                    rs_records = all_rs_variants[rs]

                    rs_without_map_weight = get_rs_without_map_weight(rs_records)
                    if not rs_without_map_weight:
                        logger.error(f"All variants for {rs} are map-weighted : \n{rs_records}")
                        continue
                    elif len(rs_without_map_weight) > 1:
                        logger.error(f"More than one variant without map-weight found for RS {rs} :"
                                     f"\n {rs_without_map_weight}")
                        continue

                    ss_records = all_ss_variants[rs]
                    if not ss_records:
                        logger.error(f"No original SS record found for RS {rs}, /n {ss_records}")
                        continue

                    rs_variant = rs_records[0]

                    if check_all_ss_has_same_info(ss_records):
                        if rs_variant['start'] == ss_records[0]['start']:
                            logger.error(f"RS's and original SS's Start matches. Nothing to do")
                            continue

                        process_all_ss_has_same_info(rs_variant, ss_records)

                    else:
                        logger.error(f"Not all original SS has same info : \n {ss_records}")
                        process_ss_diff_issue_rs_split(rs_variant, ss_records)


def process_ss_diff_issue_rs_split(rs_variant, ss_records):
    # TODO : Split RS for discordant SS variants
    pass


def process_all_ss_has_same_info(rs_variant, ss_records):
    rs_with_new_start = copy.copy(rs_variant)
    rs_with_new_start['start'] = ss_records[0]['start']
    rs_with_new_start['_id'] = get_clustered_SHA1(rs_with_new_start)

    result, variant_in_db = check_for_hash_collision(rs_with_new_start['_id'])

    if result:
        logger.warn(f"Hash collision will occur for RS {rs_variant['accession']} with RS {variant_in_db['accession']}")
        resolve_collision_and_insert_rs(rs_variant, rs_with_new_start, variant_in_db)
    else:
        logger.info(f"No hash collision for RS {rs_variant['accession']}")
        correct_rs_start_and_insert(rs_variant, rs_with_new_start)


def correct_rs_start_and_insert(rs_variant, rs_with_new_start):
    dbsnp_cve_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_CLUSTERED_VARIANT_ENTITY]

    drop_statements = []
    insert_statements = []
    drop_statements.append(pymongo.DeleteOne({'_id': rs_variant['_id']}))
    insert_statements.append(pymongo.InsertOne(rs_with_new_start))

    # result_drop = dbsnp_cve_collection.bulk_write(requests=drop_statements, ordered=False)
    # logger.info('There was %s old documents dropped' % result_drop.deleted_count)
    # result_insert = dbsnp_cve_collection.bulk_write(requests=insert_statements, ordered=False)
    # logger.info('There was %s new documents inserted' % result_insert.inserted_count)


def resolve_collision_and_insert_rs(rs_variant, rs_with_new_start, variant_in_db):
    dbsnp_cve_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_CLUSTERED_VARIANT_ENTITY]
    dbsnp_cvoe_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY]

    dbsnp_cve_drop_statements = []
    dbsnp_cve_insert_statements = []
    dbsnp_cvoe_insert_statements = []

    if rs_with_new_start['accession'] < variant_in_db['accession']:
        dbsnp_cve_drop_statements.append(pymongo.DeleteOne({'_id': rs_variant['_id']}))
        dbsnp_cve_drop_statements.append(pymongo.DeleteOne({'_id': variant_in_db['_id']}))
        dbsnp_cve_insert_statements.append(pymongo.InsertOne(rs_with_new_start))

        merge_event = create_merge_event(variant_in_db, rs_with_new_start)
        dbsnp_cvoe_insert_statements.append(pymongo.InsertOne(merge_event))

        # result_drop = dbsnp_cve_collection.bulk_write(requests=dbsnp_cve_drop_statements, ordered=False)
        # logger.info('There was %s old documents dropped' % result_drop.deleted_count)
        # result_cve_insert = dbsnp_cve_collection.bulk_write(requests=dbsnp_cve_insert_statements, ordered=False)
        # logger.info('There was %s new documents inserted' % result_cve_insert.inserted_count)
        # result_cvoe_insert = dbsnp_cvoe_collection.bulk_write(requests=dbsnp_cvoe_insert_statements, ordered=False)
        # logger.info('There was %s new documents inserted' % result_cvoe_insert.inserted_count)

        update_ss_with_new_rs(variant_in_db['accession'], rs_with_new_start['accession'])

    else:
        dbsnp_cve_drop_statements.append(pymongo.DeleteOne({'_id': rs_variant['_id']}))

        # result_drop = dbsnp_cve_collection.bulk_write(requests=dbsnp_cve_drop_statements, ordered=False)
        # logger.info('There was %s old documents dropped' % result_drop.deleted_count)

        update_ss_with_new_rs(rs_with_new_start['accession'], variant_in_db['accession'])


def create_merge_event(variant_merged, variant_retained):
    merge_event = {
        "eventType": "MERGED",
        "accession": variant_merged['accession'],
        "mergeInto": variant_retained['accession'],
        "reason": "Merged after RS discordant correction",
        "inactiveObjects": [variant_merged]
    }

    return merge_event


def update_ss_with_new_rs(old_rs, new_rs):
    dbsnp_sve_collection = mongo_source.mongo_handle[mongo_source.db_name][DBSNP_SUBMITTED_VARIANT_ENTITY]
    eva_sve_collection = mongo_source.mongo_handle[mongo_source.db_name][EVA_SUBMITTED_VARIANT_ENTITY]

    filter_query = {'rs': old_rs}
    update_value = {'$set': {'rs': new_rs}}

    # dbsnp_result_updated = dbsnp_sve_collection.update_many(filter_query, update_value)
    # logger.info('There was %s documents dbsnp updated' % dbsnp_result_updated.modified_count)
    # eva_result_updated = eva_sve_collection.update_many(filter_query, update_value)
    # logger.info('There was %s documents eva updated' % eva_result_updated.modified_count)


def get_rs_variants(mongo_source, assembly, rs_list):
    rs_filter_criteria = {'asm': assembly, 'accession': {'$in': rs_list}}
    dbsnp_rs_variants = get_variants(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY, rs_filter_criteria, 'accession')
    eva_rs_variants = get_variants(mongo_source, EVA_CLUSTERED_VARIANT_ENTITY, rs_filter_criteria, 'accession')
    all_rs_variants = merge_all_records(dbsnp_rs_variants, eva_rs_variants)
    return dbsnp_rs_variants, eva_rs_variants, all_rs_variants


def get_ss_variants(mongo_source, assembly, rs_list):
    ss_filter_criteria = {'seq': assembly, 'rs': {'$in': rs_list}}
    dbsnp_ss_variants = get_variants(mongo_source, DBSNP_SUBMITTED_VARIANT_ENTITY, ss_filter_criteria, 'rs')
    eva_ss_variants = get_variants(mongo_source, EVA_SUBMITTED_VARIANT_ENTITY, ss_filter_criteria, 'rs')
    all_ss_variants = merge_all_records(dbsnp_ss_variants, eva_ss_variants)
    return dbsnp_ss_variants, eva_ss_variants, all_ss_variants


def get_rs_events(mongo_source, rs_list):
    event_filter_criteria = {'$or': [{'accession': {'$in': rs_list}},
                                     {'mergeInto': {'$in': rs_list}},
                                     {'splitInto': {'$in': rs_list}}]}
    dbsnp_events = get_events(mongo_source, DBSNP_CLUSTERED_VARIANT_OPERATION_ENTITY, event_filter_criteria)
    eva_events = get_events(mongo_source, EVA_CLUSTERED_VARIANT_OPERATION_ENTITY, event_filter_criteria)
    all_events = merge_all_records(dbsnp_events, eva_events)
    return dbsnp_events, eva_events, all_events


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
    cursor = collection.with_options(read_concern=ReadConcern("majority")) \
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


def check_all_ss_has_same_info(ss_list):
    contig = set([ss['contig'] for ss in ss_list])
    start = set([ss['start'] for ss in ss_list])

    if len(contig) > 1 or len(start) > 1:
        return False
    else:
        return True


def check_for_hash_collision(id):
    hash_collision_filter_criteria = {'_id': id}
    collision_records = find_documents(mongo_source, DBSNP_CLUSTERED_VARIANT_ENTITY, hash_collision_filter_criteria)

    if collision_records:
        return True, collision_records[0]
    else:
        return False, None


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
                                 db_name="eva_accession_sharded")

    all_files = [os.path.join(args.discordant_rs_dir, filename) for filename in os.listdir(args.discordant_rs_dir)]
    for file in sorted(all_files, key=lambda x: os.stat(x).st_size):
        assembly = os.path.basename(file)
        fix_discordant_variants(mongo_source, assembly, os.path.join(args.discordant_rs_dir, assembly))
