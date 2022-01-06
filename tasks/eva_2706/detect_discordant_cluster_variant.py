import argparse
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


def detect_discordant_cluster_variant(mongo_source, assemblies, batch_size=1000):
    dbsnp_cve_collection = mongo_source.mongo_handle[mongo_source.db_name]["dbsnpClusteredVariantEntity"]
    dbsnp_sve_collection = mongo_source.mongo_handle[mongo_source.db_name]["dbsnpSubmittedVariantEntity"]
    filter_criteria = {}
    if assemblies:
        filter_criteria = {'asm': {"$in": assemblies}}
    cursor = dbsnp_cve_collection.with_options(read_concern=ReadConcern("majority"))\
                                 .find(filter_criteria, no_cursor_timeout=True)
    cursor.batch_size(batch_size)
    projection = {'contig': 1, 'start': 1, 'rs': 1}
    nb_clustered_variants = 0
    nb_error = 0
    for clustered_variants in grouper(batch_size, cursor):
        rsids = [clustered_variant.get('accession') for clustered_variant in clustered_variants if clustered_variant]
        nb_clustered_variants += len(rsids)

        sve_cursor = dbsnp_sve_collection.with_options(read_concern=ReadConcern("majority"))\
                                         .find({'rs': {'$in': rsids}}, projection)
        submitted_variant_position_per_rs = defaultdict(set)
        for sve in sve_cursor:
            submitted_variant_position_per_rs[sve.get('rs')].add(f"{sve.get('contig')}:{sve.get('start')}")

        nb_error += compare(clustered_variants, submitted_variant_position_per_rs)
        logger.info(f'Processed {nb_clustered_variants}: Found {nb_error} errors')


def main():
    parser = argparse.ArgumentParser(
        description='Delete declustered variants in dbsnpSubmittedVariantEntity Collection')
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
    detect_discordant_cluster_variant(mongo_source, args.assemblies, args.batch_size)


if __name__ == "__main__":
    main()
