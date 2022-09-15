from argparse import ArgumentParser
from collections import Counter, defaultdict

from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from itertools import zip_longest
from pymongo import MongoClient


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def check_submitted_variant_have_duplicates(submitted_variant_records):
    count_per_ssid = Counter()
    ssid_to_variant_list = defaultdict(list)
    duplicated_variants = []
    for variant_rec in submitted_variant_records:
        if 'allelesMatch' in variant_rec or 'mapWeight' in variant_rec:
            # Ingore the variant that have multimapped or allele mismatch flag set because they are ignored anyway
            continue
        count_per_ssid[variant_rec['accession']] += 1
        ssid_to_variant_list[variant_rec['accession']].append(variant_rec)

    for ssid in count_per_ssid:
        if count_per_ssid[ssid] > 1:
            duplicated_variants.extend(ssid_to_variant_list[ssid])
    return duplicated_variants


def shelve_submitted_variant_entities(mongo_handle, submitted_variant_accession, assembly_accession):
    output_collection = 'eva2979_dbsnpSubmittedVariantEntity'
    batch_size = 1000
    for batch_sve_accs in grouper(submitted_variant_accession, batch_size):
        query_filter = {'seq': assembly_accession, 'accession': {'$in': batch_sve_accs}}
        documents = [d for d in mongo_handle["eva_accession_sharded"]['dbsnpSubmittedVariantEntity'].find(query_filter)]
        duplicated_submitted_variant_ids = check_submitted_variant_have_duplicates(documents)
        mongo_handle["eva_accession_sharded"][output_collection].insert(duplicated_submitted_variant_ids)
        response = mongo_handle['dbsnpSubmittedVariantEntity'].delete_many({'_id': {'$in': duplicated_submitted_variant_ids}})
        assert response.deleted_count == len(duplicated_submitted_variant_ids), 'Not all variants were deleted from dbsnpSubmittedVariantEntity'


def parse_duplicate_list(duplicates_file):
    ignored_type = ('In_original_assembly,Same_variants', 'Remapped,Same_variants,Same_variants_in_source')
    duplicated_accessions = []
    assembly_accession = None
    with open(duplicates_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split('\t')
            if sp_line[2] not in ignored_type:
                duplicated_accessions.append(sp_line[0])
                if assembly_accession:
                    assert assembly_accession == sp_line[1]
                assembly_accession = sp_line[1]

    return duplicated_accessions, assembly_accession


def main():
    parser = ArgumentParser()
    parser.add_argument('--duplicates_file', required=True)
    parser.add_argument('--settings_xml_file', required=True)
    parser.add_argument('--profile', default='development')
    args = parser.parse_args()
    mongo_uri = get_mongo_uri_for_eva_profile(args.profile, args.settings_xml_file)
    with MongoClient(mongo_uri) as mongo_handle:
        duplicated_accessions, assembly_accession = parse_duplicate_list(args.duplicates_file)
        shelve_submitted_variant_entities(mongo_handle, duplicated_accessions, assembly_accession)


if __name__ == '__main__':
    main()
