import argparse
import os
from datetime import datetime
from itertools import islice

from bson import ObjectId
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()

batch_size = 10000
assembly = "GCA_000298735.1"


def find_ids_of_remapped_submitted_variants(mongo_source, output_dir):
    collections = ["dbsnpSubmittedVariantEntity", "submittedVariantEntity"]
    filter_criteria = {'remappedFrom': {'$exists': True}, 'seq': assembly}
    return extract_ids_to_file(mongo_source, collections, filter_criteria, output_dir)


def find_ids_of_remapped_clustered_variants(mongo_source, output_dir):
    collections = ["dbsnpClusteredVariantEntity", "clusteredVariantEntity"]
    filter_criteria = {'asm': assembly, 'createdDate': {'$gt': datetime.strptime("2021-11-10", '%Y-%m-%d')}}
    return extract_ids_to_file(mongo_source, collections, filter_criteria, output_dir)


def find_ids_of_submitted_variant_operations(mongo_source, output_dir):
    collections = ["submittedVariantOperationEntity"]
    filter_criteria = {'inactiveObjects.seq': assembly, 'eventType': {'$in': ["RS_MERGE_CANDIDATES", "RS_SPLIT_CANDIDATES"]}}
    '''"inactiveObjects.seq": "GCA_000298735.1", eventType: {$in: ["RS_MERGE_CANDIDATES", "RS_SPLIT_CANDIDATES"]}}'''
    return extract_ids_to_file(mongo_source, collections, filter_criteria, output_dir)


def extract_ids_to_file(mongo_source, collections, filter_criteria, output_dir):
    file_names = []
    for collection_name in collections:
        file_name = f"{output_dir}/{collection_name}_{assembly}.txt"
        file_names.append(file_name)
        logger.info(f'Searching in collections {collection_name}')
        collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
        cursor = collection \
            .with_options(read_concern=ReadConcern("majority")) \
            .find(filter_criteria, no_cursor_timeout=True)
        with open(file_name, "w") as file:
            for variant in cursor:
                file.write(str(variant['_id']) + '\n')
    return file_names


def delete_variants(mongo_source, files_with_ids):
    for file_with_ids in files_with_ids:
        collection_name = os.path.basename(file_with_ids).split('_')[0]
        collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
        logger.info(f"""variants deletion in process for file {file_with_ids}""")
        with open(file_with_ids, 'r') as file:
            id_count = 0
            deletion_count = 0
            while True:
                batch_ids = list(islice(file, batch_size))
                if not batch_ids:
                    break
                # remove trailing \n
                if 'Operation' in collection_name:
                    batch_ids = [ObjectId(i.strip()) for i in batch_ids]
                else:
                    batch_ids = [i.strip() for i in batch_ids]
                x = collection.delete_many({'_id': {'$in': batch_ids}})
                id_count += len(batch_ids)
                deletion_count += len(x.deleted_count)

            logger.info(f"""{id_count} id read and {deletion_count} documents deleted from collection {collection_name}""")


def main():
    parser = argparse.ArgumentParser(
        description='Delete document associated with Ovis aries remapped data in multiple collections', add_help=False)
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
    deletion_files = []
    deletion_files.extend(find_ids_of_remapped_submitted_variants(mongo_source, args.output_dir))
    deletion_files.extend(find_ids_of_remapped_clustered_variants(mongo_source, args.output_dir))
    deletion_files.extend(find_ids_of_submitted_variant_operations(mongo_source, args.output_dir))

    delete_variants(mongo_source, deletion_files)


if __name__ == "__main__":
    main()
