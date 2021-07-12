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


def get_SHA1(variant_rec):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    h = hashlib.sha1()
    keys = ['seq', 'study', 'contig', 'start', 'ref', 'alt']
    h.update('_'.join([str(variant_rec[key]) for key in keys]).encode())
    return h.hexdigest().upper()


def swap_with_correct_contig(mongo_source, contig_swap_list):
    sve_collection = mongo_source.mongo_handle[mongo_source.db_name]["submittedVariantEntity"]
    insert_statements = []
    drop_statements = []
    for contig in contig_swap_list:
        contig_insert_stmt_1, contig_drop_stmt_1 = get_insert_statements(sve_collection, contig['contig_1'],
                                                                         contig['contig_2'])
        contig_insert_stmt_2, contig_drop_stmt_2 = get_insert_statements(sve_collection, contig['contig_2'],
                                                                         contig['contig_1'])
        insert_statements.extend(contig_insert_stmt_1)
        insert_statements.extend(contig_insert_stmt_2)
        drop_statements.extend(contig_drop_stmt_1)
        drop_statements.extend(contig_drop_stmt_2)

    total_inserted = 0
    total_dropped = 0
    try:
        result_insert = sve_collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
            .bulk_write(requests=insert_statements, ordered=False)
        total_inserted += result_insert.inserted_count
        logger.info('%s / %s new documents inserted' % (total_inserted, 146))
        result_drop = sve_collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
            .bulk_write(requests=drop_statements, ordered=False)
        total_dropped += result_drop.deleted_count
        logger.info('%s / %s old documents dropped' % (total_dropped, 146))
    except Exception as e:
        print(traceback.format_exc())
        raise e


def get_insert_statements(sve_collection, curr_contig, swap_contig):
    filter_criteria = {'seq': 'GCA_000003025.6', 'study': 'PRJEB28579', 'contig': curr_contig}
    cursor = sve_collection.with_options(read_concern=ReadConcern("majority")) \
        .find(filter_criteria, no_cursor_timeout=True)
    insert_statements = []
    drop_statements = []
    try:
        for variant in cursor:
            original_id = get_SHA1(variant)
            assert variant['_id'] == original_id, "Original id is different from the one calculated %s != %s" % (
                variant['_id'], original_id)
            variant['contig'] = swap_contig
            variant['_id'] = get_SHA1(variant)
            insert_statements.append(pymongo.InsertOne(variant))
            drop_statements.append(pymongo.DeleteOne({'_id': original_id}))
    except Exception as e:
        print(traceback.format_exc())
        raise e
    finally:
        cursor.close()

    return insert_statements, drop_statements


def main():
    parser = argparse.ArgumentParser(description='Correct 98 contig error in study PRJEB28579', add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")

    contig_swap_list = [{'contig_1': 'AEMK02000229.1', 'contig_2': 'AEMK02000626.1'},
                        {'contig_1': 'AEMK02000417.1', 'contig_2': 'AEMK02000654.1'}]
    swap_with_correct_contig(mongo_source, contig_swap_list)


if __name__ == "__main__":
    main()