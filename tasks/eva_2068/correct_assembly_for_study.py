#!/usr/bin/env python
from argparse import ArgumentParser

from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle


def correct(mongo_user, mongo_password, mongo_host, study, reference_source, reference_dest):
    with get_mongo_connection_handle(
            username=mongo_user,
            password=mongo_password,
            host=mongo_host
    ) as accessioning_mongo_handle:
        results = accessioning_mongo_handle["eva_accession_sharded"]["submittedVariantEntity"].update_many(
            filter={'study': study, 'seq': reference_source},
            update={"$set": {'seq': reference_dest}},
            upsert=False
        )
        print(results.raw_result)
        print('There was %s documents matching and %s documents updated' % (results.matched_count, results.modified_count))


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
