#!/usr/bin/env python
import hashlib
from argparse import ArgumentParser

import pymongo
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle


def get_SHA1(variant_rec):
    """Calculate the SHA1 digest from the refee the """
    h = hashlib.sha1()
    keys = ['seq', 'study', 'contig', 'start', 'ref', 'alt']
    h.update('_'.join([str(variant_rec[key]) for key in keys]).encode())
    return h.hexdigest().upper()


def correct(mongo_user, mongo_password, mongo_host, study, reference_source, reference_dest):
    with get_mongo_connection_handle(
            username=mongo_user,
            password=mongo_password,
            host=mongo_host
    ) as accessioning_mongo_handle:
        sve_collection = accessioning_mongo_handle["eva_accession_sharded"]["submittedVariantEntity"]
        cursor = sve_collection.find({'study':study, 'seq': reference_source})
        insert_statements = []
        drop_statements = []
        record_checked = 0
        for variant in cursor:
            # Ensure that the variant we are changing has the expected SHA1
            original_id = get_SHA1(variant)
            assert variant['_id'] == original_id, "Original id is different from the one calculated %s != %s" % (variant['_id'], original_id)
            variant['seq'] = reference_dest
            variant['_id'] = get_SHA1(variant)
            insert_statements.append(pymongo.InsertOne(variant))
            drop_statements.append(pymongo.DeleteOne({'_id': original_id}))
            record_checked += 1

        print('Retrieved %s documents and checked matching Sha1 hash' % record_checked)

        result_insert = sve_collection.bulk_write(requests=insert_statements, ordered=False)
        print('There was %s new documents inserted' % result_insert.inserted_count)
        result_drop = sve_collection.bulk_write(requests=drop_statements, ordered=False)
        print('There was %s old documents dropped' % result_drop.deleted_count)


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--mongo_user', help='user to connect to mongodb', required=True)
    argparse.add_argument('--mongo_password', help='password to connect to mongodb', required=True)
    argparse.add_argument('--mongo_host', help='host to connect to mongodb', required=True)
    argparse.add_argument('--study', help='The study to correct', required=True)
    argparse.add_argument('--accession_source', help='the assembly accession of the entities that needs to be changed',
                          required=True)
    argparse.add_argument('--accession_dest', help='The assembly accession the entities will be changed to',
                          required=True)

    args = argparse.parse_args()
    correct(args.mongo_user, args.mongo_password, args.mongo_host, args.study, args.accession_source, args.accession_dest)


if __name__ == "__main__":
    main()
