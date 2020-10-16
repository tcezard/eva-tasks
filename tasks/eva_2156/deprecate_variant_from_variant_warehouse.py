#!/usr/bin/env python
from argparse import ArgumentParser

import pymongo

from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)


def deprecate(settings_xml_file, database_name, contigs=None):
    """
    Connect to mongodb and retrieve all variants that needs to be deprecated.
    Copy the variant in the operation collection and delete them from the submitted variant collections.
    """
    with pymongo.MongoClient(get_mongo_uri_for_eva_profile('production', settings_xml_file)) as mongo_handle:
        variant_collection = mongo_handle[database_name]['variants_2_0']

        cursor = variant_collection.find({'chr': {'$in': contigs}})
        drop_statements = []
        for variant in cursor:
            drop_statements.append(pymongo.DeleteOne({'_id': variant['_id']}))

    logger.info('Found %s variant to remove', len(drop_statements))

    result_drop = variant_collection.bulk_write(requests=drop_statements, ordered=False)
    logger.info('There was %s documents dropped from ' % result_drop.deleted_count)
    mongo_handle.close()


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--settings_xml_file', help='File containing the connection to the database', required=False)
    argparse.add_argument('--database_name', help='The name of the database from the variant warehouse',
                          required=True)
    argparse.add_argument('--contigs', help='The contigs to modify. they should be provided as they appeared in the record', nargs='+')
    args = argparse.parse_args()
    log_cfg.add_stdout_handler()

    deprecate(args.settings_xml_file, args.database_name, args.contigs)

    logger.info("Finished successfully.")


if __name__ == "__main__":
    main()
