import argparse
import os.path
import traceback
from collections import defaultdict
from itertools import islice

import pymongo
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern

logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)

SUBMITTED_VARIANT_ENTITY = "submittedVariantEntity"
DBSNP_SUBMITTED_VARIANT_ENTITY = "dbsnpSubmittedVariantEntity"
DBSNP_SUBMITTED_VARIANT_OPERATION_ENTITY = "dbsnpSubmittedVariantOperationEntity"


def retrieve_dbsnp_sve_with_allelesMatch_false(mongo_source, doc_file):
    retrieve_records_with_allelesMatch_false(mongo_source, DBSNP_SUBMITTED_VARIANT_ENTITY, doc_file)


def retrieve_sve_with_allelesMatch_false(mongo_source, doc_file):
    retrieve_records_with_allelesMatch_false(mongo_source, SUBMITTED_VARIANT_ENTITY, doc_file)


def retrieve_records_with_allelesMatch_false(mongo_source, collection_name, doc_file):
    with open(doc_file, 'w') as file:
        filter_criteria = {"allelesMatch": False}
        collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
        cursor = collection.with_options(read_concern=ReadConcern("majority"),
                                         read_preference=pymongo.ReadPreference.PRIMARY) \
            .find(filter_criteria, no_cursor_timeout=True)
        try:
            for record in cursor:
                ss_id = record['_id']
                ss_accession = record['accession']
                rs = record['rs'] if 'rs' in record else ''
                file.write(f'{ss_id},{ss_accession},{rs}\n')
        except Exception as e:
            logger.exception(traceback.format_exc())
            raise e
        finally:
            cursor.close()


def ss_with_same_accession_but_without_alleles_match_false(sve_doc_file, output_doc_file):
    with open(sve_doc_file, 'r') as sve_file:
        with open(output_doc_file, 'w') as output_file:
            batch_id = 1
            while True:
                logger.info(f'Processing SS should have Alleles Match False - Batch Id {batch_id}')
                batch_id = batch_id + 1
                ss_acc_list = [int(line.split(",")[1].strip()) for line in list(islice(sve_file, 100000))]
                if not ss_acc_list:
                    break
                filter_criteria = {'accession': {'$in': ss_acc_list}}
                records = find_documents(mongo_source, SUBMITTED_VARIANT_ENTITY, filter_criteria)
                dbsnp_records = find_documents(mongo_source, DBSNP_SUBMITTED_VARIANT_ENTITY, filter_criteria)
                records.extend(dbsnp_records)

                records = [record for record in records if 'allelesMatch' not in record]

                for record in records:
                    ss_id = record['_id']
                    ss_accession = record['accession']
                    output_file.write(f'{ss_id},{ss_accession}\n')


def get_rs_of_ss_with_alleles_match_false(sve_doc_file):
    rs_list_full = []
    with open(sve_doc_file, 'r') as doc_file:
        while True:
            rs_list = [line.split(",")[2].strip() for line in list(islice(doc_file, 1000000))]
            if not rs_list:
                break
            rs_list = [rs for rs in rs_list if rs != ""]
            if rs_list:
                rs_list_full.extend(rs_list)

    return rs_list_full


def ss_merged_into_ss_with_alleles_match_false(mongo_source, sve_doc_file, merged_doc_file):
    with open(sve_doc_file, 'r') as sve_file:
        with open(merged_doc_file, 'w') as merged_file:
            batch_id = 1
            while True:
                logger.info(f'Pocessing Merge Batch Id {batch_id}')
                batch_id = batch_id + 1
                ss_acc_list = [int(line.split(",")[1].strip()) for line in list(islice(sve_file, 100000))]
                if not ss_acc_list:
                    break
                filter_criteria = {'eventType': 'MERGED', 'mergeInto': {'$in': ss_acc_list}}
                records = find_documents(mongo_source, DBSNP_SUBMITTED_VARIANT_OPERATION_ENTITY, filter_criteria)
                for record in records:
                    event_id = record['_id']
                    merged_into_ss_accession = record['mergeInto']
                    merged_ss_accession = record['accession']
                    merged_ss_id = record['inactiveObjects'][0]['hashedMessage']
                    merged_ss_rs = record['inactiveObjects'][0]['rs'] if 'rs' in record['inactiveObjects'][0] else ''
                    merged_file.write(
                        f'{event_id},{merged_into_ss_accession},{merged_ss_accession},{merged_ss_id},{merged_ss_rs}\n')


def check_all_mergedInto_entities_has_alleles_match_false(merged_doc_file, output_file):
    with open(merged_doc_file, 'r') as input_file:
        with open(output_file, 'w') as output_file:
            batch_id = 1
            while True:
                logger.info(f'Processing Merge Batch Id {batch_id}')
                batch_id = batch_id + 1
                lines_list = list(islice(input_file, 100000))
                if not lines_list:
                    break

                mergeInto_acc_list = [int(line.split(",")[1].strip()) for line in lines_list]
                filter_criteria = {'accession': {'$in': mergeInto_acc_list}}
                records = find_documents(mongo_source, SUBMITTED_VARIANT_ENTITY, filter_criteria)
                dbsnp_records = find_documents(mongo_source, DBSNP_SUBMITTED_VARIANT_ENTITY, filter_criteria)
                records.extend(dbsnp_records)
                merged_acc_records = defaultdict(list)
                for record in records:
                    record_acc = record['accession']
                    merged_acc_records[record_acc].append(record)

                for line in lines_list:
                    mergee_acc = int(line.split(",")[1].strip())
                    mergee_hash = line.split(",")[3].strip()
                    found = False
                    if mergee_acc in merged_acc_records:
                        for record in merged_acc_records[mergee_acc]:
                            if record['_id'] == mergee_hash:
                                if 'allelesMatch' in record and record['allelesMatch'] == False:
                                    found = True
                                    break
                            else:
                                logger.info(f"record hash not found: {line}")

                    if not found:
                        output_file.write(f"{line}")


def ss_split_from_ss_with_alleles_match_false(mongo_source, sve_doc_file, split_doc_file):
    with open(sve_doc_file, 'r') as sve_file:
        with open(split_doc_file, 'w') as split_file:
            batch_id = 1
            while True:
                ss_acc_list = [int(line.split(",")[1].strip()) for line in list(islice(sve_file, 100000))]
                logger.info(f'Processing Split Batch Id {batch_id}')
                batch_id = batch_id + 1
                if not ss_acc_list:
                    break
                filter_criteria = {'eventType': 'SS_SPLIT', 'accession': {'$in': ss_acc_list}}
                records = find_documents(mongo_source, DBSNP_SUBMITTED_VARIANT_OPERATION_ENTITY, filter_criteria)
                for record in records:
                    event_id = record['_id']
                    split_from_ss_accession = record['accession']
                    split_into_ss_accession = record["splitInto"]
                    split_file.write(f'{event_id},{split_from_ss_accession},{split_into_ss_accession}\n')


def check_if_splitInto_ss_has_alleles_match_false(split_doc_file, split_output_file):
    with open(split_doc_file, 'r') as sve_input_file:
        with open(split_output_file, 'w') as sve_output_file:
            batch_id = 1
            while True:
                ss_acc_list = [int(line.split(",")[2].strip()) for line in list(islice(sve_input_file, 100000))]
                logger.info(f'Processing Split Batch Id {batch_id}')
                batch_id = batch_id + 1
                if not ss_acc_list:
                    break
                filter_criteria = {'accession': {'$in': ss_acc_list}}
                records = find_documents(mongo_source, SUBMITTED_VARIANT_ENTITY, filter_criteria)
                dbsnp_records = find_documents(mongo_source, DBSNP_SUBMITTED_VARIANT_ENTITY, filter_criteria)
                records.extend(dbsnp_records)

                records = [record for record in records if 'allelesMatch' in record]

                for record in records:
                    ss_id = record['_id']
                    ss_accession = record['accession']
                    sve_output_file.write(f"{ss_id},{ss_accession},{record['allelesMatch']}\n")


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count Variant stats with allelesMatch=False', add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)",
                        required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--res-dir", help="File containing discordant rs ids", required=True)
    args = parser.parse_args()

    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")

    # ##### no document with allelesMatch=false in submittedVariantEntity
    # sve_doc_file = os.path.join(args.res_dir, "sve_records_alleles_match_false")
    # if not os.path.isfile(sve_doc_file):
    #     retrieve_sve_with_allelesMatch_false(mongo_source, sve_doc_file)

    dbsnp_sve_doc_file = os.path.join(args.res_dir, "dbsnp_sve_records_alleles_match_false")
    if not os.path.isfile(dbsnp_sve_doc_file):
        retrieve_dbsnp_sve_with_allelesMatch_false(mongo_source, dbsnp_sve_doc_file)

    sve_should_have = os.path.join(args.res_dir, "ss_that_should_have_allele_match_false")
    ss_with_same_accession_but_without_alleles_match_false(dbsnp_sve_doc_file, sve_should_have)

    rs_list = get_rs_of_ss_with_alleles_match_false(dbsnp_sve_doc_file)
    if not rs_list:
        logger.info("None of the dbsnp sve retrieved has rs associated with it")

    merged_doc_file = os.path.join(args.res_dir, "ss_merged_into_ss_with_alleles_match_false")
    ss_merged_into_ss_with_alleles_match_false(mongo_source, dbsnp_sve_doc_file, merged_doc_file)

    output_file = os.path.join(args.res_dir, "not_mergedInto_having_alleles_match_false")
    check_all_mergedInto_entities_has_alleles_match_false(merged_doc_file, output_file)

    split_doc_file = os.path.join(args.res_dir, "ss_split_from_ss_with_alleles_match_false")
    ss_split_from_ss_with_alleles_match_false(mongo_source, dbsnp_sve_doc_file, split_doc_file)

    splitInto_ss_with_alleles_match_false = os.path.join(args.res_dir, "splitInto_ss_with_alleles_match_False")
    check_if_splitInto_ss_has_alleles_match_false(split_doc_file, splitInto_ss_with_alleles_match_false)

    logger.info(f"Process Finished")
