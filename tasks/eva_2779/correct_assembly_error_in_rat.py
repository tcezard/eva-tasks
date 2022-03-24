import argparse
import hashlib
import traceback

import pymongo
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo import WriteConcern
from pymongo.read_concern import ReadConcern

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


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
            result_insert = collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000))\
                                      .bulk_write(requests=insert_statements, ordered=False)
            total_inserted += result_insert.inserted_count
            result_drop = collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
                                    .bulk_write(requests=drop_statements, ordered=False)
            total_dropped += result_drop.deleted_count
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
                filter_dict.update(filter_criteria)
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
            result_update = collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
                .bulk_write(requests=update_statements, ordered=False)
            total_updated += result_update.modified_count
        logger.info(f'{total_updated} documents updated in {collection_name}')
    except Exception as e:
        print(traceback.format_exc())
        raise e
    return total_updated


def replace_variant_entities(mongo_source, batch_size):
    # submitted variants on MT chromosome
    replace_document_with_correct_information(
        mongo_source, "submittedVariantEntity", get_submitted_SHA1,
        {'seq': 'GCA_015227675.1', 'contig': 'CM026996.1'},
        {'seq': 'GCA_015227675.2', 'contig': 'AY172581.1'},
        batch_size
    )
    replace_document_with_correct_information(
        mongo_source, "dbsnpSubmittedVariantEntity", get_submitted_SHA1,
        {'seq': 'GCA_015227675.1', 'contig': 'CM026996.1'},
        {'seq': 'GCA_015227675.2', 'contig': 'AY172581.1'},
        batch_size
    )

    # submitted variants Not on MT chromosome
    replace_document_with_correct_information(
        mongo_source, "submittedVariantEntity", get_submitted_SHA1, {'seq': 'GCA_015227675.1'},
        {'seq': 'GCA_015227675.2'}, batch_size
    )
    replace_document_with_correct_information(
        mongo_source, "dbsnpSubmittedVariantEntity", get_submitted_SHA1, {'seq': 'GCA_015227675.1'},
        {'seq': 'GCA_015227675.2'}, batch_size
    )

    # Clustered variants on MT chromosome
    replace_document_with_correct_information(
        mongo_source, "clusteredVariantEntity", get_clustered_SHA1,
        {'asm': 'GCA_015227675.1', 'contig': 'CM026996.1'},
        {'asm': 'GCA_015227675.2', 'contig': 'AY172581.1'},
        batch_size
    )
    replace_document_with_correct_information(
        mongo_source, "dbsnpClusteredVariantEntity", get_clustered_SHA1,
        {'asm': 'GCA_015227675.1', 'contig': 'CM026996.1'},
        {'asm': 'GCA_015227675.2', 'contig': 'AY172581.1'},
        batch_size
    )

    # Clustered variants Not on MT chromosome
    replace_document_with_correct_information(
        mongo_source, "clusteredVariantEntity", get_clustered_SHA1, {'asm': 'GCA_015227675.1'},
        {'asm': 'GCA_015227675.2'}, batch_size
    )
    replace_document_with_correct_information(
        mongo_source, "dbsnpClusteredVariantEntity", get_clustered_SHA1, {'asm': 'GCA_015227675.1'},
        {'asm': 'GCA_015227675.2'}, batch_size
    )

    # MERGED Operations on MT chromosomes
    update_operation_entities(
        mongo_source, 'dbsnpClusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'inactiveObjects.contig': 'CM026996.1', 'eventType': 'MERGED'},
        {'inactiveObjects.asm': 'GCA_015227675.2', 'inactiveObjects.contig': 'AY172581.1', 'reason': lambda x: x.replace('GCA_015227675.1', 'GCA_015227675.2')},
        batch_size
    )
    update_operation_entities(
        mongo_source, 'clusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'inactiveObjects.contig': 'CM026996.1', 'eventType': 'MERGED'},
        {'inactiveObjects.asm': 'GCA_015227675.2', 'inactiveObjects.contig': 'AY172581.1', 'reason': lambda x: x.replace('GCA_015227675.1', 'GCA_015227675.2')},
        batch_size
    )

    # MERGED Operations Not on MT chromosomes
    update_operation_entities(
        mongo_source, 'dbsnpClusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'eventType': 'MERGED'},
        {'inactiveObjects.asm': 'GCA_015227675.2', 'reason': lambda x: x.replace('GCA_015227675.1', 'GCA_015227675.2')},
        batch_size
    )
    update_operation_entities(
        mongo_source, 'clusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'eventType': 'MERGED'},
        {'inactiveObjects.asm': 'GCA_015227675.2', 'reason': lambda x: x.replace('GCA_015227675.1', 'GCA_015227675.2')},
        batch_size
    )

    # SPLIT Operations on MT chromosomes
    update_operation_entities(
        mongo_source, 'dbsnpClusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'inactiveObjects.contig': 'CM026996.1', 'eventType': 'RS_SPLIT'},
        {'inactiveObjects.asm': 'GCA_015227675.2', 'inactiveObjects.contig': 'AY172581.1'},
        batch_size
    )
    update_operation_entities(
        mongo_source, 'clusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'inactiveObjects.contig': 'CM026996.1', 'eventType': 'RS_SPLIT'},
        {'inactiveObjects.asm': 'GCA_015227675.2', 'inactiveObjects.contig': 'AY172581.1'},
        batch_size
    )

    # SPLIT Operations Not on MT chromosomes
    update_operation_entities(
        mongo_source, 'dbsnpClusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'eventType': 'RS_SPLIT'},
        {'inactiveObjects.asm': 'GCA_015227675.2'},
        batch_size
    )
    update_operation_entities(
        mongo_source, 'clusteredVariantOperationEntity', get_clustered_SHA1,
        {'inactiveObjects.asm': 'GCA_015227675.1', 'eventType': 'RS_SPLIT'},
        {'inactiveObjects.asm': 'GCA_015227675.2'},
        batch_size
    )

    # UPDATE Operations on MT chromosome
    update_operation_entities(
        mongo_source, 'dbsnpSubmittedVariantOperationEntity', get_submitted_SHA1,
        {'inactiveObjects.seq': 'GCA_015227675.1', 'inactiveObjects.contig': 'CM026996.1', 'eventType': 'UPDATED'},
        {'inactiveObjects.seq': 'GCA_015227675.2', 'inactiveObjects.contig': 'AY172581.1'},
        batch_size
    )
    update_operation_entities(
        mongo_source, 'submittedVariantOperationEntity', get_submitted_SHA1,
        {'inactiveObjects.seq': 'GCA_015227675.1', 'inactiveObjects.contig': 'CM026996.1', 'eventType': 'UPDATED'},
        {'inactiveObjects.seq': 'GCA_015227675.2', 'inactiveObjects.contig': 'AY172581.1'},
        batch_size
    )

    # UPDATE Operations Not on MT chromosome
    update_operation_entities(
        mongo_source, 'dbsnpSubmittedVariantOperationEntity', get_submitted_SHA1,
        {'inactiveObjects.seq': 'GCA_015227675.1', 'eventType': 'UPDATED'},
        {'inactiveObjects.seq': 'GCA_015227675.2'},
        batch_size
    )
    update_operation_entities(
        mongo_source, 'submittedVariantOperationEntity', get_submitted_SHA1,
        {'inactiveObjects.seq': 'GCA_015227675.1', 'eventType': 'UPDATED'},
        {'inactiveObjects.seq': 'GCA_015227675.2'},
        batch_size
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
    parser.add_argument("--batch-size", help="number of document processed at once", required=False, default=1000)
    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")
    replace_variant_entities(mongo_source, batch_size=args.batch_size)


if __name__ == "__main__":
    main()
