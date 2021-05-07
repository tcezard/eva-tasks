import argparse
from argparse import ArgumentParser
from datetime import datetime

import psycopg2
import psycopg2.extras
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from ebi_eva_common_pyutils.pg_utils import execute_query
from retry import retry

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()
mongo_migration_count_validation_table_name = "eva_tasks.mongo4_migration_count_validation"


def create_collection_count_validation_report(mongo_source: MongoDatabase, database_list):
    count_validation_res_list = []
    report_timestamp = datetime.now()
    mongo_host = mongo_source.mongo_handle.address[0]

    for db in database_list:
        mongo_source.db_name = db
        source_collections = mongo_source.get_collection_names()

        if not source_collections:
            logger.warning(f"database {db} does not exist in mongo instances {mongo_host}")
            continue

        for coll in sorted(source_collections):
            logger.info(f"fetching count for database ({db}) - collection ({coll})")

            no_of_documents = get_documents_count_for_collection(mongo_source, db, coll)
            logger.info(f"Found {no_of_documents} documents in database ({db}) - collection ({coll})")

            count_validation_res_list.append([mongo_host, db, coll, no_of_documents, report_timestamp])

    return count_validation_res_list


@retry(logger=logger, tries=3, delay=3, backoff=2)
def get_documents_count_for_collection(mongo_server: MongoDatabase, db, coll):
    return mongo_server.mongo_handle[db][coll].count_documents({})


def create_table_for_count_validation(private_config_xml_file):
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                          user="evadev") as metadata_connection_handle:
        query_create_table_for_count_validation = "create table if not exists {0} " \
                                                  "(mongo_host text, database text, collection text, " \
                                                  "document_count bigint not null, report_time timestamp, " \
                                                  "primary key(mongo_host, database, collection, report_time))" \
            .format(mongo_migration_count_validation_table_name)

    execute_query(metadata_connection_handle, query_create_table_for_count_validation)


def insert_count_validation_result_to_db(private_config_xml_file, count_validation_res_list):
    if len(count_validation_res_list) > 0:
        with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                              user="evadev") as metadata_connection_handle:
            with metadata_connection_handle.cursor() as cursor:
                psycopg2.extras.execute_values(cursor,
                                               "INSERT INTO {0} "
                                               "(mongo_host, database, collection, document_count,report_time) "
                                               "VALUES %s".format(mongo_migration_count_validation_table_name),
                                               [tuple(x) for x in count_validation_res_list])


def get_databases_list_for_validation(file_path):
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
    parser: ArgumentParser = argparse.ArgumentParser(
        description='Given a list of databases, calculate number of documents in each collection of each of the given database and store the report',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--db-list",
                        help="Full path to the File containing list of Databases (ex: eva_hsapiens_grch37)",
                        required=True)
    parser.add_argument("--private-config-xml-file",
                        help="ex: /path/to/eva-maven-settings.xml",
                        required=True)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file)
    database_list = get_databases_list_for_validation(args.db_list)

    create_table_for_count_validation(args.private_config_xml_file)
    count_validation_result = create_collection_count_validation_report(mongo_source, database_list)
    insert_count_validation_result_to_db(args.private_config_xml_file, count_validation_result)


if __name__ == "__main__":
    main()
