import argparse

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle

from tasks.eva_2641.constants import annotation_metadata_collection_name, temp_collection_name

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def cleanup_metadata(settings_xml_file, db_name):
    with get_mongo_connection_handle('production', settings_xml_file) as mongo_conn:
        db = mongo_conn[db_name]
        metadata_collection = db[annotation_metadata_collection_name]
        temp_collection = db[temp_collection_name]

        results = [x for x in metadata_collection.find({'ct': {'$exists': True}})]
        insert_result = temp_collection.insert_many(results)
        logger.info(f'Inserted {len(insert_result.inserted_ids)} documents into {temp_collection_name}')
        delete_result = metadata_collection.delete_many(results)
        logger.info(f'Deleted {delete_result.deleted_count} documents from {annotation_metadata_collection_name}')
        metadata_collection.drop_indexes()
        logger.info(f'Dropped non-id indexes from {annotation_metadata_collection_name}')


def main():
    parser = argparse.ArgumentParser(
        description='Cleanup metadata collection by removing annotation documents to a temporary collection and '
                    'removing unnecessary indexes', add_help=False)
    parser.add_argument('--settings-xml-file', required=True)
    parser.add_argument('--db-name', required=True)

    args = parser.parse_args()
    cleanup_metadata(args.settings_xml_file, args.db_name)


if __name__ == '__main__':
    main()
