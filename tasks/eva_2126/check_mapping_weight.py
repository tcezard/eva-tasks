#!/usr/bin/env python
from argparse import ArgumentParser
from urllib.parse import quote_plus

import pymongo


def get_mongo_connection_handle(host, port=27017, username=None, password=None, authentication_database="admin", **kwargs) -> pymongo.MongoClient:
    mongo_connection_uri = "mongodb://"
    if username and password:
        mongo_connection_uri += '%s:%s@' % (username, quote_plus(password))
    mongo_connection_uri += '%s:%s/%s' % (host, port, authentication_database)
    return pymongo.MongoClient(mongo_connection_uri, **kwargs)


def check_mapping_weight(mongo_host, database_name, username, password, assembly_accession):
    """
    Connect to mongodb and retrieve all clustered variants of specific assembly that have high mapping weight (>1) to check
    if they can be found multiple times. When they can't, check that they fall in one of the following categories:
     - One or several clustered variants have been merged with another variant
     - One or several submitted variants have been declustered (because its definition in dbsnp was inconsistent) leaving
     a clustered variant without evidence
    """
    with get_mongo_connection_handle(mongo_host, username=username, password=password) as accessioning_mongo_handle:
        dbsnp_cve_collection = accessioning_mongo_handle[database_name]["dbsnpClusteredVariantEntity"]
        cursor = dbsnp_cve_collection.find({'asm': assembly_accession, 'mapWeight': {'$gt': 1}})
        count_clustered_variants = 0
        list_variants_inconsistent = []
        inconsistent_recovered_merged = 0
        inconsistent_recovered_declusted = 0
        for variant in cursor:
            count = dbsnp_cve_collection.count({'asm': assembly_accession, 'accession': variant['accession']})
            count_clustered_variants += 1
            if count < 2:
                list_variants_inconsistent.append(variant)

        dbsnp_cve_op_collection = accessioning_mongo_handle[database_name]["dbsnpClusteredVariantOperationEntity"]
        dbsnp_sve_op_collection = accessioning_mongo_handle[database_name]["dbsnpSubmittedVariantOperationEntity"]

        for variant in list_variants_inconsistent:
            merged_cve_operation = list(dbsnp_cve_op_collection.find({'accession': variant['accession'], 'eventType': 'MERGED'}))
            updated_sve_operation = list(dbsnp_sve_op_collection.find({'inactiveObjects.rs': variant['accession'], 'eventType': 'UPDATED'}))

            if merged_cve_operation:
                inconsistent_recovered_merged += 1
            elif updated_sve_operation and any([op['reason'].startswith('Declustered:') for op in updated_sve_operation]):
                inconsistent_recovered_declusted += 1
            else:
                print(
                    "rs%s Should have more than 1 variant because mapWeight = %s" % (
                    str(variant['accession']), variant['mapWeight'])
                )
                print("Found %s operation for clustered variants" % len(merged_cve_operation))
                print("Found %s operation for submitted variants" % len(updated_sve_operation))

        print("Checked %s clustered variants" % count_clustered_variants)
        print("Found %s inconsistent variants" % len(list_variants_inconsistent))
        print("Found %s inconsistent variants that are due to merged" % inconsistent_recovered_merged)
        print("Found %s inconsistent variants that are due to declustering" % inconsistent_recovered_declusted)


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--host', help='', required=True)
    argparse.add_argument('--database_name', help='', required=True)
    argparse.add_argument('--username', help='', default=None)
    argparse.add_argument('--password', help='', default=None)
    argparse.add_argument('--assembly_accession', help='', required=True)
    args = argparse.parse_args()

    check_mapping_weight(args.host, args.database_name, args.username, args.password, args.assembly_accession)


if __name__ == "__main__":
    main()
