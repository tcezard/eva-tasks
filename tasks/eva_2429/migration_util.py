import os.path

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase

logger = logging_config.get_logger(__name__)


def mongo_import_from_dir(mongo_dest_uri, mongo_dest_secrets_file, export_dir):
    mongo_import_args = {
        "mode": "upsert"
    }
    db_list = os.listdir(export_dir)

    for db in db_list:
        mongo_dest = MongoDatabase(uri=mongo_dest_uri, secrets_file=mongo_dest_secrets_file, db_name=db)
        db_dir = os.path.join(export_dir, db)
        all_coll_dir = os.listdir(db_dir)
        for coll in all_coll_dir:
            logger.info(f'Importing data for db ({db} - collection ({coll})')
            coll_dir = os.path.join(db_dir, coll)
            files_list = os.listdir(coll_dir)
            for file in files_list:
                mongo_import_args.update({"collection": coll})
                mongo_dest.import_data(os.path.join(coll_dir, file), mongo_import_args)


def write_query_to_file(query, query_file_dir, file_name):
    if query_file_dir and os.path.exists(query_file_dir):
        query_file = os.path.join(query_file_dir, file_name)
        logger.info(f"Creating query file in the location : {query_file}")
    else:
        query_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)
        logger.info(f"No query file path provided. Creating query file in the default location : {query_file}")
    try:
        with open(query_file, 'w') as open_file:
            open_file.write(query)
    except IOError:
        logger.error(f"Could not write query to the file: {query_file}")
    else:
        logger.info(f"Query written successfully to file: {query_file}")

    return query_file
