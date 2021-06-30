import argparse
import hashlib
import logging
import traceback

import pymongo
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo import WriteConcern
from pymongo.read_concern import ReadConcern


def get_SHA1(variant_rec):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    h = hashlib.sha1()
    keys = ['seq', 'study', 'contig', 'start', 'ref', 'alt']
    h.update('_'.join([str(variant_rec[key]) for key in keys]).encode())
    return h.hexdigest().upper()


def replace_with_correct_contig(mongo_source):
    correct_contig = 'AF034253.1'
    sve_collection = mongo_source.mongo_handle[mongo_source.db_name]["submittedVariantEntity"]
    filter_criteria = {'seq': 'GCA_000003025.4', 'study': 'PRJEB43246', 'contig': 'M'}
    cursor = sve_collection.with_options(read_concern=ReadConcern("majority")) \
        .find(filter_criteria, no_cursor_timeout=True)
    insert_statements = []
    drop_statements = []
    number_of_variants_to_replace = 10
    total_inserted, total_dropped = 0, 0
    try:
        for variant in cursor:
            original_id = get_SHA1(variant)
            assert variant['_id'] == original_id, "Original id is different from the one calculated %s != %s" % (
                variant['_id'], original_id)
            variant['contig'] = correct_contig
            variant['_id'] = get_SHA1(variant)
            insert_statements.append(pymongo.InsertOne(variant))
            drop_statements.append(pymongo.DeleteOne({'_id': original_id}))
        result_insert = sve_collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
            .bulk_write(requests=insert_statements, ordered=False)
        total_inserted += result_insert.inserted_count
        result_drop = sve_collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
            .bulk_write(requests=drop_statements, ordered=False)
        total_dropped += result_drop.deleted_count
        logging.info('%s / %s new documents inserted' % (total_inserted, number_of_variants_to_replace))
        logging.info('%s / %s old documents dropped' % (total_dropped, number_of_variants_to_replace))
    except Exception as e:
        print(traceback.format_exc())
        raise e
    finally:
        cursor.close()
    return total_inserted


def main():
    parser = argparse.ArgumentParser(
        description='Correct contig error in study PRJEB43246 by replacing contig M with AF034253.1', add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")
    replace_with_correct_contig(mongo_source)


if __name__ == "__main__":
    main()
