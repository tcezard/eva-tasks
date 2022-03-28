import argparse
import sys

from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern


def find_rs_entity_not_exist_in_collection(mongo_db, collection_name, rsid_list, assembly_accession):
    cve_collection = mongo_db.mongo_handle[mongo_db.db_name][collection_name]
    filter_criteria = {'asm': assembly_accession, 'accession': {'$in': rsid_list}}
    projection = {'accession': 1}
    cursor = cve_collection.with_options(read_concern=ReadConcern("majority"))\
                           .find(filter_criteria, projection, no_cursor_timeout=True)
    accession_found = [record['accession'] for record in cursor]
    remaining_rsids = set(rsid_list).difference(accession_found)
    return list(remaining_rsids)


def find_rs_entity_not_exist(mongo_db, collections, rsid_list, assembly_accession):
    remaining_rsids = find_rs_entity_not_exist_in_collection(mongo_db, collections[0], rsid_list, assembly_accession)
    if remaining_rsids:
        remaining_rsids = find_rs_entity_not_exist_in_collection(mongo_db, collections[1], remaining_rsids, assembly_accession)
    return remaining_rsids


def find_rs_references_in_ss_collection(mongo_db, collection_name, assembly_accession, batch_size=1000):
    sve_collection = mongo_db.mongo_handle[mongo_db.db_name][collection_name]
    filter_criteria = {'seq': assembly_accession, 'rs': {'$exists': True}}
    projection = {'rs': 1}
    cursor = sve_collection.with_options(read_concern=ReadConcern("majority"))\
                           .find(filter_criteria, projection, no_cursor_timeout=True)
    rs_list = []
    for record in cursor:
        rs_list.append(record['rs'])
        if len(rs_list) == batch_size:
            yield rs_list
            rs_list = []
    yield rs_list


def check_rs_for_assembly(mongo_db, assembly_accession, batch_size):
    nb_processed = 0
    for rs_list in find_rs_references_in_ss_collection(mongo_db, 'submittedVariantEntity', assembly_accession, batch_size):
        if rs_list:
            nb_processed += len(rs_list)
            print(f'Processes {nb_processed} EVA submitted variants', file=sys.stderr)
            rs_with_no_entity = find_rs_entity_not_exist(
                mongo_db, ['clusteredVariantEntity', 'dbsnpClusteredVariantEntity'], rs_list, assembly_accession
            )
            if rs_with_no_entity:
                for rs in rs_with_no_entity:
                    print(f'Found a EVA submitted variant entity referencing rs {rs} but no clustered variant entity was found for it.')
    nb_processed = 0
    for rs_list in find_rs_references_in_ss_collection(mongo_db, 'dbsnpSubmittedVariantEntity', assembly_accession, batch_size):
        if rs_list:
            nb_processed += len(rs_list)
            print(f'Processes {nb_processed} DBSNP submitted variants', file=sys.stderr)
            rs_with_no_entity = find_rs_entity_not_exist(
                mongo_db, ['dbsnpClusteredVariantEntity', 'clusteredVariantEntity'], rs_list, assembly_accession
            )
            if rs_with_no_entity:
                for rs in rs_with_no_entity:
                    print(f'Found a DBSNP submitted variant entity referencing rs {rs} but no clustered variant entity was found for it.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check if the RS entities referenced in SS entities all exist in '
                                                 'the database')
    parser.add_argument("--assembly_accession",
                        help="The assembly accession used to retrieve submitted and clustered variant entities", required=True)
    parser.add_argument("--mongo-db-uri",
                        help="Mongo Database URI (ex: mongodb://user:@mongos-db-host:27017/admin)", required=True)
    parser.add_argument("--mongo-db-secrets-file",
                        help="Full path to the Mongo Database secrets file (ex: /path/to/mongo/db/secret)",
                        required=True)
    parser.add_argument("--batch_size", help="Number of SS queried in one batch", default=2000, required=False)

    args = parser.parse_args()
    mongo_db = MongoDatabase(uri=args.mongo_db_uri, secrets_file=args.mongo_db_secrets_file,
                             db_name="eva_accession_sharded")
    check_rs_for_assembly(mongo_db, args.assembly_accession, args.batch_size)
