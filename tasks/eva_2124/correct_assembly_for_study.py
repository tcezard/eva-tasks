#!/usr/bin/env python
import hashlib
from argparse import ArgumentParser
import pymongo
from urllib.parse import quote_plus
from tasks.eva_2124.load_synonyms import load_synonyms_for_assembly


def get_SHA1(variant_rec):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    h = hashlib.sha1()
    keys = ['seq', 'study', 'contig', 'start', 'ref', 'alt']
    h.update('_'.join([str(variant_rec[key]) for key in keys]).encode())
    return h.hexdigest().upper()


def get_genbank(synonym_dictionaries, contig):
    """
    returns a tuple (genbank, was_already_genbank) or raises an exception if the contig was not found
    """
    by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc = synonym_dictionaries
    if contig in by_name:
        return by_name[contig]['genbank'], False

    if contig in by_assigned_molecule:
        return by_assigned_molecule[contig]['genbank'], False

    if contig in by_ucsc:
        return by_ucsc[contig]['genbank'], False

    if contig in by_refseq and by_refseq[contig]['is_genbank_refseq_identical']:
        return by_refseq[contig]['genbank'], False

    if contig in by_genbank:
        return contig, True

    raise Exception('could not find synonym for contig {}'.format(contig))


def get_mongo_connection_handle_url(host, port=27017, username=None, password=None, authentication_database="admin") -> pymongo.MongoClient:
    mongo_connection_uri = "mongodb://"
    if username and password:
        mongo_connection_uri += '%s:%s@' % (username, quote_plus(password))
    mongo_connection_uri += '%s:%s/%s' % (host, port, authentication_database)
    return pymongo.MongoClient(mongo_connection_uri)


def correct(mongo_user, mongo_password, mongo_host, studies, assembly_accession, mongo_database, assembly_report=None,
            chunk_size=1000):
    """
    Connect to mongodb and retrieve all variants the should be updated, check their key and update them in bulk.
    """
    with get_mongo_connection_handle_url(
            username=mongo_user,
            password=mongo_password,
            host=mongo_host
    ) as accessioning_mongo_handle:
        sve_collection = accessioning_mongo_handle[mongo_database]["submittedVariantEntity"]
        synonym_dictionaries = load_synonyms_for_assembly(assembly_accession, assembly_report)
        assert_all_contigs_can_be_replaced(sve_collection, synonym_dictionaries, studies, assembly_accession)
        return do_updates(sve_collection, synonym_dictionaries, studies, assembly_accession, chunk_size)


def assert_all_contigs_can_be_replaced(sve_collection, synonym_dictionaries, studies, assembly_accession):
    cursor = sve_collection.aggregate([{'$match': {'seq': assembly_accession, 'study': {'$in': studies}}},
                                       {'$group': {'_id': '$contig', 'count': {'$sum': 1}}}])
    already_genbank_contigs = 0
    already_genbank_variants = 0
    replaceable_contigs = 0
    replaceable_variants = 0
    unreplaceable_variants = 0
    unreplaceable_contigs = set()
    for contig in cursor:
        try:
            genbank, was_already_genbank = get_genbank(synonym_dictionaries, contig['_id'])
            if was_already_genbank:
                already_genbank_contigs += 1
                already_genbank_variants += contig['count']
            else:
                replaceable_contigs += 1
                replaceable_variants += contig['count']
        except Exception:
            unreplaceable_variants += contig['count']
            unreplaceable_contigs.add(contig['_id'])

    if len(unreplaceable_contigs) > 0:
        raise Exception(
            'Aborting replacement (no changes were done).'
            ' With the provided assembly_report, the next {} contigs (present in {} variants) can not be replaced: {}'
            .format(len(unreplaceable_contigs), unreplaceable_variants, unreplaceable_contigs))
    elif replaceable_contigs == 0:
        raise Exception('All the variants were already genbank ({} variants, {} contigs). '
                        'Are you sure the assembly ({}) and studies ({}) are correct?'
                        .format(already_genbank_variants, already_genbank_contigs, assembly_accession, studies))
    else:
        print('Check ok. {} contigs of {} variants will be changed to genbank, '
              'and {} contigs of {} variants are already genbank contigs.'
              .format(replaceable_contigs, replaceable_variants, already_genbank_contigs, already_genbank_variants))


def do_updates(sve_collection, synonym_dictionaries, studies, assembly_accession, chunk_size):
    cursor = sve_collection.find({'study': {'$in': studies}, 'seq': assembly_accession})
    insert_statements = []
    drop_statements = []
    record_checked = 0
    already_genbanks = 0
    total_inserted = 0
    total_dropped = 0
    print("Performing updates...")
    for variant in cursor:
        # Ensure that the variant we are changing has the expected SHA1
        original_id = get_SHA1(variant)
        assert variant['_id'] == original_id, "Original id is different from the one calculated %s != %s" % (
            variant['_id'], original_id)
        genbank, was_already_genbank = get_genbank(synonym_dictionaries, variant['contig'])
        if was_already_genbank:
            already_genbanks += 1
        else:
            variant['contig'] = genbank
            variant['_id'] = get_SHA1(variant)
            insert_statements.append(pymongo.InsertOne(variant))
            drop_statements.append(pymongo.DeleteOne({'_id': original_id}))
        record_checked += 1
        if len(insert_statements) >= chunk_size:
            total_inserted, total_dropped = execute_bulk(drop_statements, insert_statements, sve_collection,
                                                         total_dropped, total_inserted)
    if len(insert_statements) > 0:
        total_inserted, total_dropped = execute_bulk(drop_statements, insert_statements, sve_collection,
                                                     total_dropped, total_inserted)
    print('Retrieved %s documents and checked matching Sha1 hash' % record_checked)
    print('{} of those documents had already a genbank contig. If the projects were all affected, '
          'that number should be 0, but even if it is not, there is nothing else to fix'.format(already_genbanks))
    print('There was %s new documents inserted' % total_inserted)
    print('There was %s old documents dropped' % total_dropped)
    return total_inserted


def execute_bulk(drop_statements, insert_statements, sve_collection, total_dropped, total_inserted):
    result_insert = sve_collection.bulk_write(requests=insert_statements, ordered=False)
    total_inserted += result_insert.inserted_count
    result_drop = sve_collection.bulk_write(requests=drop_statements, ordered=False)
    total_dropped += result_drop.deleted_count
    insert_statements.clear()
    drop_statements.clear()
    return total_inserted, total_dropped


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--mongo_user', help='user to connect to mongodb', required=False)
    argparse.add_argument('--mongo_password', help='password to connect to mongodb', required=False)
    argparse.add_argument('--mongo_host', help='host to connect to mongodb', required=True)
    argparse.add_argument('--mongo_database', help='The DB where the submittedVariantEntity is stored',
                          required=False, default='eva_accession_sharded')
    argparse.add_argument('--studies', help='The studies in the assembly to correct', required=True, nargs='+')
    argparse.add_argument('--assembly', help='The assembly accession of the entities that needs to be changed',
                          required=True)
    argparse.add_argument('--assembly-report', help='Use this assembly report instead of downloading the standard one',
                          required=False)

    args = argparse.parse_args()
    correct(args.mongo_user, args.mongo_password, args.mongo_host, args.studies, args.assembly, args.mongo_database,
            args.assembly_report)


if __name__ == "__main__":
    main()
