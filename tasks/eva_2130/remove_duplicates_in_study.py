#!/usr/bin/env python
from argparse import ArgumentParser
import pymongo
from urllib.parse import quote_plus


def same_variant_except_contig(variant1, variant2):
    list_key = ['tax', 'study', 'start', 'ref', 'alt', 'version', 'seq']
    return all([variant1[key] == variant2[key] for key in list_key]) and variant1['contig'] != variant2['contig']


def find_duplicates_and_remove_them(mongo_user, mongo_password, mongo_host, mongo_database, assembly_accession,
                                    contig_list, study_list, dry_run):
    duplicates_to_remove_commands = []
    with get_mongo_connection_handle_url(
            username=mongo_user,
            password=mongo_password,
            host=mongo_host
    ) as accessioning_mongo_handle:
        sve_collection = accessioning_mongo_handle[mongo_database]["submittedVariantEntity"]
        try:
            cursor = sve_collection.find({
                "seq": assembly_accession, "contig": {"$in": contig_list}, "study": {"$in": study_list}
            }, no_cursor_timeout=True)
            for record in cursor:
                variants_with_accessions = list(accessioning_mongo_handle[mongo_database]["submittedVariantEntity"].find({'accession': record['accession']}))
                if len(variants_with_accessions) != 1:
                    print('Found %s duplicates for accession %s' % (len(variants_with_accessions), record['accession']) )
                for variant in variants_with_accessions:
                    if same_variant_except_contig(variant, record):
                        duplicates_to_remove_commands.append(
                            pymongo.DeleteOne({'contig': variant['contig'], 'accession': variant['accession']})
                        )
        except Exception as e:
            raise e
        finally:
            cursor.close()

        if dry_run:
            print("Will remove %s variants " % len(duplicates_to_remove_commands))
        else:
            sve_collection.bulk_write(requests=duplicates_to_remove_commands, ordered=False)
            print("Remove %s variants " % len(duplicates_to_remove_commands))


def get_mongo_connection_handle_url(host, port=27017, username=None, password=None, authentication_database="admin") -> pymongo.MongoClient:
    mongo_connection_uri = "mongodb://"
    if username and password:
        mongo_connection_uri += '%s:%s@' % (username, quote_plus(password))
    mongo_connection_uri += '%s:%s/%s' % (host, port, authentication_database)
    return pymongo.MongoClient(mongo_connection_uri)


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--mongo_user', help='user to connect to mongodb', required=False)
    argparse.add_argument('--mongo_password', help='password to connect to mongodb', required=False)
    argparse.add_argument('--mongo_host', help='host to connect to mongodb', required=True)
    argparse.add_argument('--mongo_database', help='The DB where the submittedVariantEntity is stored',
                          required=False, default='eva_accession_sharded')
    argparse.add_argument('--assembly', help='The assembly accession of the entities that needs to be changed',
                          required=True)
    argparse.add_argument('--contigs', help='The contigs in the assembly to correct', required=True, nargs='+')
    argparse.add_argument('--studies', help='The studies in the assembly to correct', required=True, nargs='+')
    argparse.add_argument('--dry_run', help='Check that the variant contig names can be replaced using the assembly report',
                          default=False, action='store_true')

    args = argparse.parse_args()
    find_duplicates_and_remove_them(args.mongo_user, args.mongo_password, args.mongo_host, args.mongo_database, args.assembly,
                                    args.contigs, args.studies, args.dry_run)


if __name__ == "__main__":
    main()
