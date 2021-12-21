import argparse
import os
from itertools import islice

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from pymongo.read_concern import ReadConcern

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()

batch_size = 1000


def find_ids_of_declustered_variants(mongo_source, output_dir):
    affected_assemblies = [
        "GCA_001522545.2", "GCA_900700415.1", "GCA_003254395.2", "GCA_003957565.2",
        "GCA_000188115.3", "GCA_000219495.2", "GCA_000512255.2", "GCA_000001515.5",
        "GCA_000003195.3", "GCA_011100615.1", "GCA_002880775.3", "GCA_014441545.1",
        "GCA_000003025.6", "GCA_000751015.1", "GCA_002114115.1", "GCA_000181335.4",
        "GCA_000002315.5", "GCA_000146605.4", "GCA_000372685.2", "GCA_015227675.1",
        "GCA_002863925.1", "GCA_000298735.1", "GCA_001433935.1", "GCA_000002035.4",
        "GCA_001858045.3", "GCA_000001215.4", "GCA_000004515.4", "GCA_902167145.1",
        "GCA_000001635.9", "GCA_002263795.2", "GCA_001704415.1", "GCA_000002775.3",
        "GCA_000224145.1", "GCA_003339765.3", "GCA_008746955.1"
    ]

    print(f"""Number of Affected Assemblies :  {len(affected_assemblies)}""")
    dbsnp_sve_collection = mongo_source.mongo_handle[mongo_source.db_name]["dbsnpSubmittedVariantEntity"]

    for assembly in affected_assemblies:
        print('Running for assembly: ' + assembly)

        filter_criteria = {'remappedFrom': {'$exists': True}, 'rs': {'$exists': False}, 'seq': assembly}
        cursor = dbsnp_sve_collection.with_options(read_concern=ReadConcern("majority")) \
            .find(filter_criteria, no_cursor_timeout=True)

        f = open(f"""{output_dir}/{assembly}.txt""", "a")
        for variant in cursor:
            f.write(variant['_id'] + '\n')
        f.close()


def delete_variants(mongo_source, output_dir):
    dbsnp_sve_collection = mongo_source.mongo_handle[mongo_source.db_name]["dbsnpSubmittedVariantEntity"]
    for variant_file in os.listdir(output_dir):
        print(f"""variants deletion in process for assembly  {variant_file}""")
        with open(os.path.join(output_dir, variant_file), 'r') as file:
            variant_count = 0
            while True:
                batch_ids = list(islice(file, batch_size))
                if not batch_ids:
                    break
                # remove trailing \n
                batch_ids = [i.strip() for i in batch_ids]
                variant_count = variant_count + len(batch_ids)
                #dbsnp_sve_collection.deleteMany({'_id': {'$in': batch_ids}})
            print(f"""variants deleted for assembly {variant_file} : {variant_count}""")



def main():
    parser = argparse.ArgumentParser(
        description='Delete declustered variants in dbsnpSubmittedVariantEntity Collection', add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--output-dir", help="Top-level directory where all files reside (ex: /path/to/files)",
                        required=True)
    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")
    # this method finds the variants to delete and stored their ids in a file named by assembly
    find_ids_of_declustered_variants(mongo_source, args.output_dir)
    # this method reads the variant ids store by previous method in batches and deletes them
    #delete_variants(mongo_source, args.output_dir)


if __name__ == "__main__":
    main()
