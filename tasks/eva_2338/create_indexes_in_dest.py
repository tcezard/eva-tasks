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
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from ebi_eva_common_pyutils.logger import logging_config

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def create_indexes(mongo_source: MongoDatabase, mongo_dest: MongoDatabase):
    logger.info(f"Creating indexes in the target database {mongo_dest.uri_with_db_name}....")
    mongo_dest.create_index_on_collections(mongo_source.get_indexes())


def main():
    parser = argparse.ArgumentParser(description='Dump data from a given MongoDB source',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--mongo-dest-uri",
                        help="Mongo Destination URI (ex: mongodb://user:@mongos-dest-host:27017/admin)",
                        required=True)
    parser.add_argument("--mongo-dest-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--db-name", help="Database to migrate (ex: eva_hsapiens_grch37)", required=True)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()
    create_indexes(
        MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file, db_name=args.db_name),
        MongoDatabase(uri=args.mongo_dest_uri, secrets_file=args.mongo_dest_secrets_file, db_name=args.db_name))


if __name__ == "__main__":
    main()
