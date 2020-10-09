#!/usr/bin/env python
from argparse import ArgumentParser
from copy import deepcopy
from datetime import datetime

import pymongo

from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)


def inactive_object(variant):
    inactive_variant = deepcopy(variant)
    inactive_variant['hashedMessage'] = inactive_variant.pop('_id')
    return {
        "eventType": "DEPRECATED",
        "accession": variant['accession'],
        "reason": "Removed: The contig on which the variant was mapped is not accessioned by INSDC.",
        "inactiveObjects": [inactive_variant],
        "createdDate": datetime.now()  # "2018-07-30T13:21:59.333Z"
    }


def deprecate(settings_xml_file, study, assembly_accession, contigs=None):
    """
    Connect to mongodb and retrieve all variants the should be updated, check their key and update them in bulk.
    """
    with pymongo.MongoClient(get_mongo_uri_for_eva_profile('production', settings_xml_file)) as accessioning_mongo_handle:
        sve_collection = accessioning_mongo_handle['eva_accession_sharded']["submittedVariantEntity"]
        deprecated_sve_collection = accessioning_mongo_handle['eva_accession_sharded']["submittedVariantOperationEntity"]
        cursor = sve_collection.find({'seq': assembly_accession, 'study': study, 'contig': {'$in': contigs}})
        insert_statements = []
        drop_statements = []
        for variant in cursor:
            insert_statements.append(pymongo.InsertOne(inactive_object(variant)))
            drop_statements.append(pymongo.DeleteOne({'_id': variant['_id']}))

    # There should only be 458 variant to deprecate
    assert len(insert_statements) == 458
    assert len(drop_statements) == 458

    logger.info('Found %s variant to deprecate', len(insert_statements))

    result_insert = deprecated_sve_collection.bulk_write(requests=insert_statements, ordered=False)
    result_drop = sve_collection.bulk_write(requests=drop_statements, ordered=False)
    logger.info('There was %s new documents inserted in inactive entities' % result_insert.inserted_count)
    logger.info('There was %s old documents dropped from ' % result_drop.deleted_count)
    accessioning_mongo_handle.close()


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--settings_xml_file', help='File containing the connection to the database', required=False)
    argparse.add_argument('--study', help='The study in the assembly to correct', required=True)
    argparse.add_argument('--assembly', help='The assembly accession of the entities that needs to be changed',
                          required=True)
    argparse.add_argument('--contigs', help='The contigs to modify. they should be provided as they appeared in the record', nargs='+')
    args = argparse.parse_args()
    log_cfg.add_stdout_handler()

    deprecate(args.settings_xml_file, args.study, args.assembly, args.contigs)

    logger.info("Finished successfully.")


if __name__ == "__main__":
    main()
