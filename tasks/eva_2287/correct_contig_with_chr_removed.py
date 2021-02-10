import hashlib
import argparse
import pymongo
import traceback
import logging
from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile


def get_SHA1(variant_rec):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    h = hashlib.sha1()
    keys = ['seq', 'study', 'contig', 'start', 'ref', 'alt']
    h.update('_'.join([str(variant_rec[key]) for key in keys]).encode())
    return h.hexdigest().upper()


def correct(private_config_xml_file, profile='production', mongo_database='eva_accession_sharded'):
    with pymongo.MongoClient(get_mongo_uri_for_eva_profile(profile, private_config_xml_file)) as mongo_handle:
        sve_collection = mongo_handle[mongo_database]["submittedVariantEntity"]
        filter_criteria = {'seq': 'GCA_002742125.1', 'study': 'PRJEB42582'}
        cursor = sve_collection.find(filter_criteria)
        insert_statements = []
        drop_statements = []
        number_of_variants_to_replace = 10
        total_inserted, total_dropped = 0, 0
        try:
            for variant in cursor:
                original_id = get_SHA1(variant)
                assert variant['_id'] == original_id, "Original id is different from the one calculated %s != %s" % (
                    variant['_id'], original_id)
                variant['contig'] = 'CM008482.1'
                variant['_id'] = get_SHA1(variant)
                insert_statements.append(pymongo.InsertOne(variant))
                drop_statements.append(pymongo.DeleteOne({'_id': original_id}))
            result_insert = sve_collection.bulk_write(requests=insert_statements, ordered=False)
            total_inserted += result_insert.inserted_count
            result_drop = sve_collection.bulk_write(requests=drop_statements, ordered=False)
            total_dropped += result_drop.deleted_count
            logging.info('%s / %s new documents inserted' % (total_inserted, number_of_variants_to_replace))
            logging.info('%s / %s old documents dropped' % (total_dropped, number_of_variants_to_replace))
        except Exception as e:
            print(traceback.format_exc())
            raise e
        finally:
            cursor.close()
        return total_inserted


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get stats from variant warehouse', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    args = parser.parse_args()
    correct(args.private_config_xml_file)
