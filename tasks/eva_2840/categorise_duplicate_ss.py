import argparse
from collections import defaultdict
from itertools import zip_longest

from ebi_eva_common_pyutils.mongodb import MongoDatabase


def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def categorise_batch_duplicate_ss(mongo_db, ssids, assembly_accession):
    sve_collection = mongo_db.mongo_handle[mongo_db.db_name]['dbsnpSubmittedVariantEntity']
    cursor = sve_collection.find({'seq': assembly_accession, 'accession': {'$in': ssids}})
    ssid_to_type_set = defaultdict(list)
    for variant_rec in cursor:
        reasons = set()
        if 'allelesMatch' in variant_rec:
            reasons.add('Mismatching allele')
        if 'mapWeight' in variant_rec:
            reasons.add('Multi mapped')
        if 'remappedFrom' in variant_rec:
            reasons.add('Remapped')
        if not reasons:
            reasons.add('other')
        ssid_to_type_set[variant_rec['accession']].append(','.join(sorted(reasons)))

    return ssid_to_type_set


def categorise_all_rs(mongo_db, duplicate_ss_file, output_file, batch_size):
    all_ss_accessions = []
    assemblies = set()
    with open(duplicate_ss_file) as open_file:
        for line in open_file:
            count, assembly, ssid = line.strip().split()
            all_ss_accessions.append(ssid)
            assemblies.add(assembly)

    nb_processed = 0
    assert len(assemblies) == 1, 'Only one assembly per file is expected'
    assembly = assemblies.pop()
    with open(output_file, 'w') as open_file:
        for ssids in grouper(all_ss_accessions, batch_size):
            ssid_to_types = categorise_batch_duplicate_ss(mongo_db, ssids, assembly)
            for ssid in ssid_to_types:
                print(f"{ssid}\t{assembly}\t{','.join(ssid_to_types[ssid])}", file=open_file)
            nb_processed += len(ssids)
            print(f"{nb_processed}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check if the RS entities referenced in SS entities all exist in '
                                                 'the database')
    parser.add_argument("--duplicate_ss_file",
                        help="Input file containing the SSids that have duplicate entries in one assembly", required=True)
    parser.add_argument("--output_file",
                        help="Output file containing the duplicate SSids annotated with the type", required=True)
    parser.add_argument("--mongo-db-uri",
                        help="Mongo Database URI (ex: mongodb://user:@mongos-db-host:27017/admin)", required=True)
    parser.add_argument("--mongo-db-secrets-file",
                        help="Full path to the Mongo Database secrets file (ex: /path/to/mongo/db/secret)",
                        required=True)

    args = parser.parse_args()
    mongo_db = MongoDatabase(uri=args.mongo_db_uri, secrets_file=args.mongo_db_secrets_file,
                             db_name="eva_accession_sharded")
    categorise_all_rs(mongo_db, args.missing_rs_file, args.output_file)
