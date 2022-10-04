from argparse import ArgumentParser
from collections import Counter, defaultdict

from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
from itertools import zip_longest
from pymongo import MongoClient, WriteConcern, ReadPreference
from pymongo.read_concern import ReadConcern

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


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
        else:
            logger.warning(f'Found only {count_per_ssid.get(ssid, 0)} variant for ss{ssid}')
    return duplicated_variants


def shelve_submitted_variant_entities(mongo_handle, submitted_variant_accession, assembly_accession):
    output_collection = 'eva2979_dbsnpSubmittedVariantEntity'
    batch_size = 1000
    for batch_sve_accs in grouper(submitted_variant_accession, batch_size):
        # remove None
        batch_sve_accs = [acc for acc in batch_sve_accs if acc]
        query_filter = {'seq': assembly_accession, 'accession': {'$in': batch_sve_accs}}
        documents = [d for d in mongo_handle["eva_accession_sharded"]['dbsnpSubmittedVariantEntity'].
                     with_options(read_concern=ReadConcern("majority"), read_preference=ReadPreference.PRIMARY).
                     find(query_filter)]
        logger.info(f'Found {len(documents)} documents for {len(batch_sve_accs)} accessions')
        duplicated_submitted_variant_ids = check_submitted_variant_have_duplicates(documents)
        if duplicated_submitted_variant_ids:
            mongo_handle["eva_accession_sharded"][output_collection].\
                with_options(write_concern=WriteConcern("majority")).\
                insert_many(duplicated_submitted_variant_ids)
            # response = mongo_handle['dbsnpSubmittedVariantEntity'].\
            #     with_options(write_concern=WriteConcern("majority")).\
            #     delete_many({'_id': {'$in': duplicated_submitted_variant_ids}})
            # assert response.deleted_count == len(duplicated_submitted_variant_ids), 'Not all variants were deleted from dbsnpSubmittedVariantEntity'


def parse_duplicate_list(duplicates_file):
    ignored_type = ('In_original_assembly,Same_variants', 'Remapped,Same_variants,Same_variants_in_source')
    duplicated_accessions = []
    assembly_accession = None
    with open(duplicates_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            if sp_line[2] not in ignored_type:
                duplicated_accessions.append(int(sp_line[0]))
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
        if duplicated_accessions:
            logger.info(f'Will attempt to shelf {len(duplicated_accessions)} for assembly {assembly_accession}')
            shelve_submitted_variant_entities(mongo_handle, duplicated_accessions, assembly_accession)
        else:
            logger.info(f'No duplicated variants for assembly {assembly_accession}')


if __name__ == '__main__':
    main()
