#!/usr/bin/env python
import datetime
from argparse import ArgumentParser

import pymongo


def datetime_converter(o):
    if isinstance(o, datetime.datetime):
        return o.__str__()


def get_mongo_connection_handle(host, port=27017, username=None, password=None, authentication_database="admin", **kwargs) -> pymongo.MongoClient:
    return pymongo.MongoClient(
        host=host,
        port=port,
        username=username,
        password=password,
        authSource=authentication_database,
        **kwargs
    )


def check_mapping_weight(mongo_host, database_name, username, password, reference_accession):
    """
    Connect to mongodb and retrieve all variants the should be updated, Check their key and update them in bulk.
    """
    with get_mongo_connection_handle(mongo_host) as accessioning_mongo_handle:
        dbsnp_cve_collection = accessioning_mongo_handle[database_name]["dbsnpClusteredVariantEntity"]
        cursor = dbsnp_cve_collection.find({'mapWeight': {'$gt': 1}, 'asm': reference_accession})
        count_variants = 0
        list_variant_erroneous = []
        erroneous_recovered1 = 0
        erroneous_recovered2 = 0
        for variant in cursor:
            count = dbsnp_cve_collection.count({'accession': variant['accession']})
            count_variants += 1
            if count < 2:
                list_variant_erroneous.append(variant)

        dbsnp_cve_op_collection = accessioning_mongo_handle[database_name]["dbsnpClusteredVariantOperationEntity"]
        dbsnp_sve_op_collection = accessioning_mongo_handle[database_name]["dbsnpSubmittedVariantOperationEntity"]

        for variant in list_variant_erroneous:
            operation1 = list(dbsnp_cve_op_collection.find({'accession': variant['accession']}))
            operation2 = list(dbsnp_sve_op_collection.find({'inactiveObjects.rs': variant['accession']}))

            if operation1 and operation1[0]['eventType'] == 'MERGED':
                erroneous_recovered1 += 1
            elif operation2 and operation2[0]['eventType'] == 'UPDATED':
                erroneous_recovered2 += 1
            else:
                print(
                    "rs%s Should have more than 1 variant because mapWeight = %s" % (
                    str(variant['accession']), variant['mapWeight'])
                )
                print("Found %s operation for clustered variants" % len(operation1))
                print("Found %s operation for submitted variants" % len(operation2))

        print("Checked %s variants" % count_variants)
        print("Found %s erroneous variants" % len(list_variant_erroneous))
        print("Found %s erroneous variants that are due to merged" % erroneous_recovered1)
        print("Found %s erroneous variants that are due to declustering" % erroneous_recovered2)


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--host', help='', required=True)
    argparse.add_argument('--database_name', help='', required=True)
    argparse.add_argument('--username', help='', default=None)
    argparse.add_argument('--password', help='', default=None)
    argparse.add_argument('--reference_accession', help='', required=True)
    args = argparse.parse_args()

    check_mapping_weight(args.host, args.database_name, args.username, args.password, args.reference_accession)


if __name__ == "__main__":
    main()
