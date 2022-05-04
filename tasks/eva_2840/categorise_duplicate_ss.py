import argparse
from collections import defaultdict
from itertools import zip_longest

from ebi_eva_common_pyutils.mongodb import MongoDatabase


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def categorise_batch_duplicate_ss(mongo_db, ssids, assembly_accession):
    sve_collection = mongo_db.mongo_handle[mongo_db.db_name]['dbsnpSubmittedVariantEntity']
    cursor = sve_collection.find({'seq': assembly_accession, 'accession': {'$in': ssids}})
    ssid_to_type_set = defaultdict(list)
    ssid_to_positions = defaultdict(set)
    ssid_to_changes = defaultdict(set)
    # Check each variant independently
    for variant_rec in cursor:
        reasons = set()
        position = f"{variant_rec['contig']}:{variant_rec['start']}"
        change = f"{variant_rec['ref']}-{variant_rec['alt']}"

        if 'allelesMatch' in variant_rec:
            reasons.add('Mismatching_allele')
        if 'mapWeight' in variant_rec:
            reasons.add('Multi_mapped')
        if not reasons:
            ssid_to_positions[variant_rec['accession']].add(position)
            ssid_to_changes[variant_rec['accession']].add(change)

            if 'remappedFrom' in variant_rec:
                reasons.add('Remapped')
            if not reasons:
                reasons.add('In_original_assembly')
        ssid_to_type_set[variant_rec['accession']].append(','.join(sorted(reasons)))
    cursor.close()
    # Check variants per ssids
    for accession in ssid_to_positions:
        if len(ssid_to_positions[accession]) > 1:
            ssid_to_type_set[accession].append('Multi_position_ssid')
        elif len(ssid_to_changes[accession]) > 1:
            ssid_to_type_set[accession].append('Multi_allele_ssid')
        else:
            ssid_to_type_set[accession].append('Same_variants')

    # check variants in different assembly
    ssids = list(ssid_to_positions)
    ssid_to_positions = defaultdict(set)
    ssid_to_changes = defaultdict(set)
    cursor = sve_collection.find({'seq': {'$ne': assembly_accession}, 'accession': {'$in': ssids}})
    for variant_rec in cursor:
        position = f"{variant_rec['contig']}:{variant_rec['start']}"
        change = f"{variant_rec['ref']}-{variant_rec['alt']}"
        ssid_to_positions[variant_rec['accession']].add(position)
        ssid_to_changes[variant_rec['accession']].add(change)
    cursor.close()

    for accession in ssid_to_positions:
        if len(ssid_to_positions[accession]) > 1:
            ssid_to_type_set[accession].append('Multi_position_in_source')
        elif len(ssid_to_changes[accession]) > 1:
            ssid_to_type_set[accession].append('Multi_allele_in_source')
        else:
            ssid_to_type_set[accession].append('Same_variants_in_source')
    return ssid_to_type_set


def categorise_all_ss(mongo_db, duplicate_ss_file, output_file, batch_size):
    all_ss_accessions = []
    assemblies = set()
    with open(duplicate_ss_file) as open_file:
        for line in open_file:
            if line.strip():
                count, assembly, ssid = line.strip().split()
                all_ss_accessions.append(int(ssid))
                assemblies.add(assembly)

    nb_processed = 0
    assert len(assemblies) == 1, 'Only one assembly per file is expected'
    assembly = assemblies.pop()
    print(f'{len(all_ss_accessions)} ssids to process')
    with open(output_file, 'w') as open_file:
        for ssids in grouper(all_ss_accessions, batch_size):
            ssids = [ssid for ssid in ssids if ssid is not None]
            ssid_to_types = categorise_batch_duplicate_ss(mongo_db, ssids, assembly)
            for ssid in ssid_to_types:
                print(f"{ssid}\t{assembly}\t{','.join(list(dict.fromkeys(ssid_to_types[ssid])))}", file=open_file)
            nb_processed += len(ssid_to_types)
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
    parser.add_argument("--batch_size", help="Number of ssids search at once", default=1000, required=False)

    args = parser.parse_args()
    mongo_db = MongoDatabase(uri=args.mongo_db_uri, secrets_file=args.mongo_db_secrets_file,
                             db_name="eva_accession_sharded")
    categorise_all_ss(mongo_db, args.duplicate_ss_file, args.output_file, args.batch_size)
    mongo_db.mongo_handle.close()
