import argparse
import os
from datetime import datetime
from itertools import islice

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()

batch_size = 10000

"""

Remove 123,270,212 + 5377 remapped submitted variants
Remove 25,555,262 remapped clustered variants

14,203,143


Submitted variants Operations
split and merge candidates
db.submittedVariantOperationEntity.find({"inactiveObjects.seq": "GCA_000298735.1", eventType: {$in: ["RS_MERGE_CANDIDATES", "RS_SPLIT_CANDIDATES"]}}).count()
14206877

db.submittedVariantEntity.find({"seq": "GCA_000298735.1", remappedFrom: {$exists:true}}).count()
db.dbsnpSubmittedVariantEntity.find({"seq": "GCA_000298735.1", remappedFrom: {$exists:true}}).count()
"""


def find_ids_of_remapped_submitted_variants(mongo_source, output_dir):
    assembly = "GCA_008746955.1"
    collections = ["dbsnpSubmittedVariantEntity", "submittedVariantEntity"]
    filter_criteria = {'remappedFrom': {'$exists': True}, 'seq': assembly}
    return extract_ids_to_file(mongo_source, collections, filter_criteria, assembly, output_dir)


def find_ids_of_remapped_clustered_variants(mongo_source, output_dir):
    assembly = "GCA_008746955.1"
    collections = ["dbsnpClusteredVariantEntity", "clusteredVariantEntity"]
    filter_criteria = {'asm': assembly, 'createdDate': {'$gt': datetime.strptime("2021-11-10", '%Y-%m-%d')}}
    return extract_ids_to_file(mongo_source, collections, filter_criteria, assembly, output_dir)


def find_ids_of_submitted_variant_operations(mongo_source, output_dir):
    assembly = "GCA_008746955.1"
    collections = ["submittedVariantOperationEntity"]
    filter_criteria = {'asm': assembly, 'eventType': {'$in': ["RS_MERGE_CANDIDATES", "RS_SPLIT_CANDIDATES"]}}
    return extract_ids_to_file(mongo_source, collections, filter_criteria, assembly, output_dir)


def extract_ids_to_file(mongo_source, collections, filter_criteria, assembly, output_dir):
    file_names = []
    for collection_name in collections:
        file_name = f"{output_dir}/{collection_name}_{assembly}.txt"
        file_names.append(file_name)
        logger.info(f'Searching in collections {collection_name}')
        cve_collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
        cursor = cve_collection \
            .with_options(read_concern=ReadConcern("majority")) \
            .find(filter_criteria, no_cursor_timeout=True)
        with open(file_name, "w") as file:
            for variant in cursor:
                file.write(variant['_id'] + '\n')
    return file_names


def delete_variants(mongo_source, output_dir):
    dbsnp_sve_collection = mongo_source.mongo_handle[mongo_source.db_name]["dbsnpSubmittedVariantEntity"]
    for variant_file in os.listdir(output_dir):
        logger.info(f"""variants deletion in process for assembly  {variant_file}""")
        with open(os.path.join(output_dir, variant_file), 'r') as file:
            variant_count = 0
            while True:
                batch_ids = list(islice(file, batch_size))
                if not batch_ids:
                    break
                # remove trailing \n
                batch_ids = [i.strip() for i in batch_ids]
                variant_count = variant_count + len(batch_ids)
                x = dbsnp_sve_collection.delete_many({'_id': {'$in': batch_ids}})
                logger.info(f"""variants deleted in the batch : {x.deleted_count}""")
            logger.info(f"""variants deleted for assembly {variant_file} : {variant_count}""")


def main():
    parser = argparse.ArgumentParser(
        description='Delete declustered variants in dbsnpSubmittedVariantEntity Collection', add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--output-dir", help="Top-level directory where all files reside (ex: /path/to/files)",
                        required=True)
    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")

    remapped_submitted_variants_file_names = find_ids_of_remapped_submitted_variants(mongo_source, args.output_dir)
    remapped_clustered_variants_file_names = find_ids_of_remapped_clustered_variants(mongo_source, args.output_dir)
    submitted_variant_operations_file_names = find_ids_of_submitted_variant_operations(mongo_source, args.output_dir)


if __name__ == "__main__":
    main()
