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

# Keys are collection names, Values are a tuple of (shard keys, uniqueness flag)
collections_shard_key_map = {"variants_2_0": (["chr", "start"], False), "files_2_0": (["sid", "fid", "fname"], True)}


def prepare_dest_db(mongo_source_db: MongoDatabase, mongo_dest_db: MongoDatabase):
    logger.info("Dropping target database if it already exists...")
    mongo_dest_db.drop()
    logger.info("Enabling sharding in the target database...")
    mongo_dest_db.enable_sharding()
    logger.info("Sharding collections in the target database...")
    mongo_dest_db.shard_collections(collections_shard_key_map,
                                    collections_to_shard=mongo_source_db.get_collection_names())


def main():
    parser = argparse.ArgumentParser(description='Prepare target database before loading',
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
    mongo_source_db = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                    db_name=args.db_name)
    mongo_dest_db = MongoDatabase(uri=args.mongo_dest_uri, secrets_file=args.mongo_dest_secrets_file,
                                  db_name=args.db_name)
    prepare_dest_db(mongo_source_db, mongo_dest_db)


if __name__ == "__main__":
    main()
