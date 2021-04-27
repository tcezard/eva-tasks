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


def create_collection_count_validation_report(mongo_source: MongoDatabase, mongo_dest: MongoDatabase, database_list):
    count_validation_res_list = []

    for db in database_list:
        mongo_source.db_name = db
        mongo_dest.db_name = db
        source_collections = mongo_source.get_collection_names()
        dest_collections = mongo_dest.get_collection_names()

        if not source_collections and not dest_collections:
            logger.info(f"database {db} does not exit in either of mongo instances")
            continue

        all_collections = sorted(set(source_collections).union(set(dest_collections)))

        for coll in all_collections:
            logger.info(f"fetching count for database {db} - collection {coll}")

            no_of_documents_in_src = get_documents_count_for_collection(mongo_source, db, coll)
            no_of_documents_in_dest = get_documents_count_for_collection(mongo_dest, db, coll)

            if (no_of_documents_in_src == 0) and (coll not in source_collections):
                logger.info(f"collection {coll} does not exist in database {db} in mongo source")
                no_of_documents_in_src = -1

            if (no_of_documents_in_dest == 0) and (coll not in dest_collections):
                logger.info(f"collection {coll} does not exist in database {db} in dest source")
                no_of_documents_in_dest = -1

            if no_of_documents_in_src != no_of_documents_in_dest and \
                    no_of_documents_in_src != -1 and no_of_documents_in_dest != -1:
                logger.info(f"no of documents in collection {coll} of database {db} does not match")

            count_validation_res_list.append(
                (db, coll, no_of_documents_in_src, no_of_documents_in_dest, datetime.now()))

    return count_validation_res_list


@retry(logger=logger, tries=3, delay=2)
def get_documents_count_for_collection(mongo_server: MongoDatabase, db, coll):
    return mongo_server.mongo_handle[db][coll].count_documents({})


def create_table_for_count_validation(private_config_xml_file):
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                          user="evadev") as metadata_connection_handle:
        query_create_table_for_count_validation = "create table if not exists {0} " \
                                                  "(database text, collection text, source_document_count bigint not null, " \
                                                  "destination_document_count bigint not null, report_time timestamp, " \
                                                  "primary key(database, collection, report_time))" \
            .format(mongo_migration_count_validation_table_name)

    execute_query(metadata_connection_handle, query_create_table_for_count_validation)


def insert_count_validation_result_to_db(private_config_xml_file, count_validation_res_list):
    if len(count_validation_res_list) > 0:
        with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                              user="evadev") as metadata_connection_handle:
            with metadata_connection_handle.cursor() as cursor:
                psycopg2.extras.execute_values(cursor,
                                               "INSERT INTO {0} "
                                               "(database, collection, source_document_count, destination_document_count,report_time) "
                                               "VALUES %s".format(mongo_migration_count_validation_table_name),
                                               count_validation_res_list)


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
    parser: ArgumentParser = argparse.ArgumentParser(
        description='Given two different mongo db instances and a list of databases, compare number of documents '
                    'in each of the corresponding collections for all the given dbs in two instances',
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
                        help="Full path to the Mongo Destination secrets file (ex: /path/to/mongo/dest/secret)",
                        required=True)
    parser.add_argument("--db-list",
                        help="Full path to the File containing list of Databases (ex: eva_hsapiens_grch37)",
                        required=True)
    parser.add_argument("--private-config-xml-file",
                        help="ex: /path/to/eva-maven-settings.xml",
                        required=True)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    source_mongo = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file)
    dest_mongo = MongoDatabase(uri=args.mongo_dest_uri, secrets_file=args.mongo_dest_secrets_file)
    database_list = get_databases_list_for_export(args.db_list)

    create_table_for_count_validation(args.private_config_xml_file)
    count_validation_result = create_collection_count_validation_report(source_mongo, dest_mongo, database_list)
    insert_count_validation_result_to_db(args.private_config_xml_file, count_validation_result)


if __name__ == "__main__":
    main()
