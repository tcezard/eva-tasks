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
import os
import sys
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from ebi_eva_common_pyutils.logger import logging_config

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def restore_data_to_dest(mongo_dest: MongoDatabase, top_level_dump_dir):
    try:
        dump_dir = os.path.join(top_level_dump_dir, mongo_dest.db_name)
        logger.info(f"Loading data in target database from source dump {dump_dir}...")
        # noIndexRestore - Do not restore indexes because MongoDB 3.2 does not have index compatibility with MongoDB 4.0
        mongo_dest.restore_data(dump_dir=dump_dir,
                                mongorestore_args={"noIndexRestore": "",
                                                   "numParallelCollections": 4,
                                                   "numInsertionWorkersPerCollection": 4})
    except Exception as ex:
        logger.error(f"Error while restoring data to the destination database!\n{ex.__str__()}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description='Restore data to a given MongoDB destination',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--mongo-dest-uri",
                        help="Mongo Destination URI (ex: mongodb://user:@mongos-dest-host:27017/admin)",
                        required=True)
    parser.add_argument("--mongo-dest-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--db-name", help="Database to migrate (ex: eva_hsapiens_grch37)", required=True)
    parser.add_argument("--dump-dir", help="Top-level directory where all dumps reside (ex: /path/to/dumps)",
                        required=True)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()
    restore_data_to_dest(MongoDatabase(args.mongo_dest_uri, args.mongo_dest_secrets_file, args.db_name),
                         top_level_dump_dir=args.dump_dir)


if __name__ == "__main__":
    main()
