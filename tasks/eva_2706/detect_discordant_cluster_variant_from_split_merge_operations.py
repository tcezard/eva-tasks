import argparse
import hashlib
from collections import defaultdict
from itertools import zip_longest

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def compare(clustered_variants, submitted_variant_position_per_rs):
    error = 0
    for clustered_variant in clustered_variants:
        if not clustered_variant:
            continue
        pos = f'{clustered_variant.get("contig")}:{clustered_variant.get("start")}'
        if clustered_variant['accession'] not in submitted_variant_position_per_rs:
            logger.error(f'No submitted variant found for rs{clustered_variant["accession"]}')
            error += 1
            continue

        positions = submitted_variant_position_per_rs[clustered_variant['accession']]
        if len(positions) > 1:
            logger.error(f'{len(positions)} positions found in submitted variants for '
                         f'rs{clustered_variant["accession"]}')
            error += 1
            continue

        if pos not in positions:
            logger.error(f'cluster position ({pos}) different from submitted position ({positions.pop()}) '
                         f'for rs{clustered_variant["accession"]}')
            error += 1
            continue
    return error


NO_SEQUENCE_ALTERATION = 'NO_SEQUENCE_ALTERATION'
SNV = 'SNV'
MNV = 'MNV'
DEL = 'DEL'
INS = 'INS'


def classify_variant(reference, alternate):
    reference = reference.strip().upper()
    alternate = alternate.strip().upper()

    if reference == "NOVARIATION" or alternate == reference:
        return NO_SEQUENCE_ALTERATION
    else:
        is_ref_alpha = reference.isalpha()
        is_alt_alpha = alternate.isalpha()

        if is_ref_alpha and is_alt_alpha:
            if len(reference) == len(alternate):
                return SNV if len(reference) == 1 else MNV
        if is_ref_alpha and not alternate:
            return DEL
        if is_alt_alpha and not reference:
            return INS


def submitted_variant_to_clustered_variant_hash(submitted_variant):
    """Calculate the SHA1 digest from the clustered variant based on the submitted variant information"""
    h = hashlib.sha1()
    h.update('_'.join([
        submitted_variant.get('seq'),
        submitted_variant.get('contig'),
        str(submitted_variant.get('start')),
        classify_variant(submitted_variant.get('ref'), submitted_variant.get('alt')),
    ]).encode())
    return h.hexdigest().upper()


def detect_discordant_cluster_variant_from_split_merge_operations(mongo_source, assemblies, batch_size=1000):
    sveo_collection = mongo_source.mongo_handle[mongo_source.db_name]["submittedVariantOperationEntity"]
    dbsnp_cve_collection = mongo_source.mongo_handle[mongo_source.db_name]["dbsnpClusteredVariantEntity"]
    dbsnp_sve_collection = mongo_source.mongo_handle[mongo_source.db_name]["dbsnpSubmittedVariantEntity"]
    sveo_filter_criteria = {'eventType': {'$in': ['RS_MERGE_CANDIDATES', 'RS_SPLIT_CANDIDATES']}}
    if assemblies:
        sveo_filter_criteria['inactiveObjects.seq'] = {"$in": assemblies}
    cursor = sveo_collection.with_options(read_concern=ReadConcern("majority"))\
                            .find(sveo_filter_criteria)
    cursor.batch_size(batch_size)
    projection = {'contig': 1, 'start': 1, 'rs': 1}
    nb_clustered_variants = 0
    nb_error = 0
    for submitted_variant_operations in grouper(batch_size, cursor):
        clustered_variant_hashes = set((
            submitted_variant_to_clustered_variant_hash(submitted_variant)
            for submitted_variant_operation in submitted_variant_operations
            for submitted_variant in submitted_variant_operation.get('inactiveObjects') if submitted_variant_operation
        ))
        clustered_variants = dbsnp_cve_collection.find({'id': {'$in': clustered_variant_hashes}})
        rsids = [clustered_variant.get('accession') for clustered_variant in clustered_variants if clustered_variant]
        nb_clustered_variants += len(rsids)
        sve_filtering = {'seq': {"$in": assemblies}, 'rs': {'$in': rsids}}
        sve_cursor = dbsnp_sve_collection.with_options(read_concern=ReadConcern("majority")) \
                                         .find(sve_filtering, projection)
        submitted_variant_position_per_rs = defaultdict(set)
        for sve in sve_cursor:
            submitted_variant_position_per_rs[sve.get('rs')].add(f"{sve.get('contig')}:{sve.get('start')}")

        nb_error += compare(clustered_variants, submitted_variant_position_per_rs)
        logger.info(f'Processed {nb_clustered_variants}: Found {nb_error} errors')


def main():
    parser = argparse.ArgumentParser(
        description='Detect clustered variant that have position discordant with their submitted variants and are '
                    'involve in a merge or split event.')
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--assemblies", nargs='+', help="The list of assembly to check", default=[])
    parser.add_argument("--batch_size", default=1000, help="The number of variant to retrieve pr batch")
    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")
    detect_discordant_cluster_variant_from_split_merge_operations(mongo_source, args.assemblies, args.batch_size)


if __name__ == "__main__":
    main()
