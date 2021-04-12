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

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def archive_data_from_source(mongo_source: MongoDatabase, top_level_archive_dir):
    try:
        logger.info("Running mongodump from source...")

        # Force table scan is performant for many workloads avoids cursor timeout issues
        # See https://jira.mongodb.org/browse/TOOLS-845?focusedCommentId=988298&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-988298
        mongo_source.archive_data(archive_dir=os.path.join(top_level_archive_dir, mongo_source.db_name),
                                  archive_name=mongo_source.db_name,
                                  mongodump_args={"gzip": "", "forceTableScan": "", "numParallelCollections": "1"})
    except Exception as ex:
        logger.error(f"Error while dumping data from source!\n{ex.__str__()}")
        sys.exit(1)


def get_databases_list_for_export(file_path):
    database_list = []
    try:
        with open(file_path) as file_object:
            lines = file_object.readlines()
    except FileNotFoundError:
        logger.error('Could not find file with database list to export. Please check file path and name.')
    else:
        for line in lines:
            database_list.append(line.strip())

    return database_list


def main():
    parser = argparse.ArgumentParser(description='Archive data from a given MongoDB source',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--db-names-list-file",
                        help="Full path to the File containing list of Databases to migrate (ex: eva_hsapiens_grch37)",
                        required=True)
    parser.add_argument("--archive-dir", help="Top-level directory where all archives reside (ex: /path/to/archives)",
                        required=True)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    databases_list = get_databases_list_for_export(args.db_names_list_file)

    for db in databases_list:
        archive_data_from_source(MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                               db_name=db), top_level_archive_dir=args.archive_dir)


if __name__ == "__main__":
    main()
