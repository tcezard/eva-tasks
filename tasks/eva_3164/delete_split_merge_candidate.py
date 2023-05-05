#!/usr/bin/env python3

# Copyright 2021 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import pymongo
import traceback

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo import WriteConcern
from pymongo.read_concern import ReadConcern


logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


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


def delete_split_and_merge_candidates(mongo_source, assembly_accession):
    sveo_collection = mongo_source.mongo_handle[mongo_source.db_name]["submittedVariantOperationEntity"]
    filter_criteria = {"eventType": {"$in": ["RS_MERGE_CANDIDATES", "RS_SPLIT_CANDIDATES"]}, "inactiveObjects.seq": assembly_accession}
    nb_documents, total_dropped = 0, 0
    for variant_batch in find_documents_in_batch(mongo_source, "submittedVariantOperationEntity", filter_criteria):
        drop_statements = [pymongo.DeleteOne({'_id': operation['_id']}) for operation in variant_batch]
        result_drop = sveo_collection.with_options(write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
            .bulk_write(requests=drop_statements, ordered=False)
        total_dropped += result_drop.deleted_count
        nb_documents += len(variant_batch)
        logger.info(f'Deleted {result_drop.deleted_count}')
    logger.info(f'{total_dropped} documents dropped out of {nb_documents}')


def main():
    parser = argparse.ArgumentParser(
        description='Remove SPLIT and MERGE CANDIDATE for the provided assembly', add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--assembly-accession", help="Genbank assembly accession (ex: GCA_000003055.5)", required=True)

    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")
    delete_split_and_merge_candidates(mongo_source, args.assembly_accession)


if __name__ == "__main__":
    main()