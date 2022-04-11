import argparse
import hashlib
import logging
import traceback

import pymongo
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo import WriteConcern
from pymongo.errors import BulkWriteError
from pymongo.read_concern import ReadConcern

logger = logging_config.get_logger(__name__)


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


def find_documents_in_batch(mongo_source, collection_name, filter_criteria, batch_size=1000):
    collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
    cursor = collection.with_options(read_concern=ReadConcern("majority")) \
                       .find(filter_criteria, no_cursor_timeout=True)
    records = []
    try:
        for result in cursor:
            records.append(result)
            if len(records) == batch_size:
                yield records
                records = []
        yield records
    except Exception as e:
        logger.exception(traceback.format_exc())
        raise e
    finally:
        cursor.close()


def replace_document_with_correct_information(mongo_source, collection_name, id_creation_func, filter_criteria,
                                              correction_map, batch_size):
    collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]

    total_inserted, total_dropped = 0, 0
    try:
        for batch_of_variant in find_documents_in_batch(mongo_source, collection_name, filter_criteria, batch_size):
            insert_statements = []
            drop_statements = []
            for variant in batch_of_variant:
                original_id = id_creation_func(variant)
                assert variant['_id'] == original_id, "Original id is different from the one calculated %s != %s" % (
                    variant['_id'], original_id)
                for key in correction_map:
                    variant[key] = correction_map[key]
                variant['_id'] = id_creation_func(variant)
                insert_statements.append(pymongo.InsertOne(variant))
                drop_statements.append(pymongo.DeleteOne({'_id': original_id}))
            if insert_statements and drop_statements:
                try:
                    result_insert = collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000))\
                                              .bulk_write(requests=insert_statements, ordered=False)
                    total_inserted += result_insert.inserted_count
                    logger.debug(f'{result_insert.inserted_count} new documents inserted in {collection_name}')
                except BulkWriteError as bulk_error:
                    error_code_names = set([error.get('codeName') for error in bulk_error.details.get('writeErrors')])
                    if len(error_code_names) == 1 and error_code_names.pop() == 'DuplicateKey':
                        # This error occurs because we were able to create the entry in a previous run but not able to
                        # remove the original variant yet
                        total_inserted += bulk_error.details.get('nInserted')
                        logger.debug(f"Duplicate key error found while inserting but still inserted "
                                     f"{bulk_error.details.get('nInserted')} new documents inserted "
                                     f"in {collection_name}")
                    else:
                        raise bulk_error

                result_drop = collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
                                        .bulk_write(requests=drop_statements, ordered=False)
                total_dropped += result_drop.deleted_count
                logger.debug(f'{result_drop.deleted_count} old documents dropped in {collection_name}')
            else:
                logger.warning(f'{len(insert_statements)} insert statements and {len(drop_statements)} drop statements '
                               f'created. Skipping.')

        logger.info(f'{total_inserted} new documents inserted in {collection_name}')
        logger.info(f'{total_dropped} old documents dropped in {collection_name}')
    except Exception as e:
        print(traceback.format_exc())
        raise e
    return total_inserted


def update_operation_entities(mongo_source, collection_name, id_creation_func, filter_criteria, correction_map, batch_size):
    collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
    total_updated = 0
    try:
        for batch_of_variant in find_documents_in_batch(mongo_source, collection_name, filter_criteria, batch_size):
            update_statements = []
            for variant in batch_of_variant:
                filter_dict = {'_id': variant['_id']}
                for key in correction_map:
                    if callable(correction_map[key]):
                        variant[key] = correction_map[key](variant[key])
                    elif 'inactiveObjects' in key:
                        inactive_objects = variant['inactiveObjects']
                        prop = key.split('.')[-1]
                        for inactive in inactive_objects:
                            inactive[prop] = correction_map[key]
                    else:
                        variant[key] = correction_map[key]
                for inactive in variant['inactiveObjects']:
                    inactive['hashedMessage'] = id_creation_func(inactive)
                variant.pop('_id')
                update_statements.append(pymongo.ReplaceOne(filter_dict, variant))
            if update_statements:
                result_update = collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
                    .bulk_write(requests=update_statements, ordered=False)
                total_updated += result_update.modified_count
                logger.debug(f'{result_update.modified_count} documents updated in {collection_name}')
            else:
                logger.warning(f'{len(update_statements)} update statements created. Skipping.')
        logger.info(f'{total_updated} documents updated in {collection_name}')
    except Exception as e:
        print(traceback.format_exc())
        raise e
    return total_updated


def replace_variant_entities(mongo_source, batch_size):
    source_asm = 'GCA_015227675.1'
    target_asm = 'GCA_015227675.2'
    source_mt = 'CM026996.1'
    target_mt = 'AY172581.1'
    filter_sve_with_MT = {'seq': source_asm, 'contig': source_mt}
    filter_sve_without_MT = {'seq': source_asm}
    change_sve_with_MT = {'seq': target_asm, 'contig': target_mt}
    change_sve_without_MT = {'seq': target_asm}
    submitted_variant_collections = ["submittedVariantEntity", "dbsnpSubmittedVariantEntity"]

    # submitted variants on MT chromosome
    logger.info(f'Change submitted variants on MT chromosome')
    for collection in submitted_variant_collections:
        replace_document_with_correct_information(
            mongo_source, collection, get_submitted_SHA1,
            filter_sve_with_MT, change_sve_with_MT, batch_size
        )

    # submitted variants Not on MT chromosome
    logger.info(f'Change submitted variants Not on MT chromosome')
    for collection in submitted_variant_collections:
        replace_document_with_correct_information(
            mongo_source, collection, get_submitted_SHA1,
            filter_sve_without_MT, change_sve_without_MT, batch_size
        )

    filter_cve_with_MT = {'asm': source_asm, 'contig': source_mt}
    filter_cve_without_MT = {'asm': source_asm}
    change_cve_with_MT = {'asm': target_asm, 'contig': target_mt}
    change_cve_without_MT = {'asm': target_asm}
    clustered_variant_collections = ["clusteredVariantEntity", "dbsnpClusteredVariantEntity"]
    # Clustered variants on MT chromosome
    logger.info(f'Change Clustered variants on MT chromosome')
    for collection in clustered_variant_collections:
        replace_document_with_correct_information(
            mongo_source, collection, get_clustered_SHA1,
            filter_cve_with_MT, change_cve_with_MT, batch_size
        )

    # Clustered variants Not on MT chromosome
    logger.info(f'Change Clustered variants Not on MT chromosome')
    for collection in clustered_variant_collections:
        replace_document_with_correct_information(
            mongo_source, collection, get_clustered_SHA1,
            filter_cve_without_MT, change_cve_without_MT, batch_size
        )

    filter_cvoe_with_merge_MT = {'inactiveObjects.asm': source_asm, 'inactiveObjects.contig': source_mt, 'eventType': 'MERGED'}
    filter_cvoe_with_merge_without_MT = {'inactiveObjects.asm': source_asm, 'eventType': 'MERGED'}
    change_cvoe_with_merge_MT = {'inactiveObjects.asm': target_asm, 'inactiveObjects.contig': target_mt, 'reason': lambda x: x.replace(source_asm, target_asm)}
    change_cvoe_with_merge_without_MT = {'inactiveObjects.asm': target_asm, 'reason': lambda x: x.replace(source_asm, target_asm)}
    clustered_variant_operation_collections = ['clusteredVariantOperationEntity', 'dbsnpClusteredVariantOperationEntity']

    # MERGED Operations on MT chromosomes
    logger.info(f'Change MERGED Operations on MT chromosomes')
    for collection in clustered_variant_operation_collections:
        update_operation_entities(
            mongo_source, collection, get_clustered_SHA1,
            filter_cvoe_with_merge_MT, change_cvoe_with_merge_MT, batch_size
        )

    # MERGED Operations Not on MT chromosomes
    logger.info(f'Change MERGED Operations Not on MT chromosomes')
    for collection in clustered_variant_operation_collections:
        update_operation_entities(
            mongo_source, collection, get_clustered_SHA1,
            filter_cvoe_with_merge_without_MT, change_cvoe_with_merge_without_MT, batch_size
        )

    filter_cvoe_with_MT = {'inactiveObjects.asm': source_asm, 'inactiveObjects.contig': source_mt}
    filter_cvoe_without_MT = {'inactiveObjects.asm': source_asm}
    change_cvoe_with_MT = {'inactiveObjects.asm': target_asm, 'inactiveObjects.contig': target_mt}
    change_cvoe_without_MT = {'inactiveObjects.asm': target_asm}
    # SPLIT Operations on MT chromosomes
    logger.info(f'Change SPLIT Operations on MT chromosomes')
    for collection in clustered_variant_operation_collections:
        update_operation_entities(
            mongo_source, collection, get_clustered_SHA1,
            filter_cvoe_with_MT, change_cvoe_with_MT, batch_size
        )

    # SPLIT Operations Not on MT chromosomes
    logger.info(f'Change SPLIT Operations Not on MT chromosomes')
    for collection in clustered_variant_operation_collections:
        update_operation_entities(
            mongo_source, collection, get_clustered_SHA1,
            filter_cvoe_without_MT, change_cvoe_without_MT, batch_size
        )

    filter_svoe_with_merge_MT = {'inactiveObjects.seq': source_asm, 'inactiveObjects.contig': source_mt, 'eventType': 'UPDATED'}
    filter_svoe_with_merge_without_MT = {'inactiveObjects.seq': source_asm, 'eventType': 'UPDATED'}
    change_svoe_with_merge_MT = {'inactiveObjects.seq': target_asm, 'inactiveObjects.contig': target_mt}
    change_svoe_with_merge_without_MT = {'inactiveObjects.seq': target_asm}
    submitted_variant_operation_collections = ['submittedVariantOperationEntity', 'dbsnpSubmittedVariantOperationEntity']
    # UPDATED Operations on MT chromosome
    logger.info(f'Change UPDATED Operations on MT chromosome')
    for collection in submitted_variant_operation_collections:
        update_operation_entities(
            mongo_source, collection, get_submitted_SHA1,
            filter_svoe_with_merge_MT, change_svoe_with_merge_MT, batch_size
        )

    # UPDATED Operations Not on MT chromosome
    logger.info(f'Change UPDATED Operations Not on MT chromosome')
    for collection in submitted_variant_operation_collections:
        update_operation_entities(
            mongo_source, collection, get_submitted_SHA1,
            filter_svoe_with_merge_without_MT, change_svoe_with_merge_without_MT, batch_size
        )


def main():
    parser = argparse.ArgumentParser(
        description='Correct assembly error in assembly GCA_015227675.1 by replacing it GCA_015227675.2',
        add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--batch-size", help="number of document processed at once", required=False, type=int, default=1000)
    parser.add_argument("--debug", help="Set the script to output debug message", default=False, action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging_config.add_stdout_handler(logging.DEBUG)
    else:
        logging_config.add_stdout_handler()

    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")
    replace_variant_entities(mongo_source, batch_size=int(args.batch_size))
    del mongo_source


if __name__ == "__main__":
    main()
