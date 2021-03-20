import hashlib
import argparse
import pymongo
import psycopg2
import requests
import traceback
import logging
from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile


def get_genbank_contig_alias(contig):
    url = 'http://ves-hx-de.ebi.ac.uk:8080/eva/webservices/contig-alias/v1/chromosomes/refseq/' + contig
    print(url)
    response = requests.get(url)
    genbank = response.json()['_embedded']['chromosomeEntities'][0]['genbank']
    return genbank


def get_SHA1(variant_id):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    return hashlib.sha1(variant_id.encode()).hexdigest().upper()


def get_ids(assembly, variant_query_result):
    """Get the id fields as a pair (id_variant_warehouse, sve_hash, cve_hash)"""
    # TODO: Support variants with multiple files
    study = variant_query_result['files'][0]['sid']
    genbank_chr = get_genbank_contig_alias(variant_query_result['chr'])
    start = variant_query_result['start']
    ref = variant_query_result['ref']
    alt = variant_query_result['alt']
    id_variant_warehouse = variant_query_result['_id']
    sve_hash = get_SHA1(f"{assembly}_{study}_{genbank_chr}_{start}_{ref}_{alt}")
    return id_variant_warehouse, sve_hash


def get_db_name_and_assembly_accession(private_config_xml_file):
    # metadata_handle = psycopg2.connect(get_pg_metadata_uri_for_eva_profile(
    #     "development", private_config_xml_file), user="evadev")
    return 'eva_fcatus_90', 'GCA_000181335.4'
    # return 'eva_fcatus_80', 'GCA_000181335.3'
    # return 'eva_cporcellus_30', 'GCA_000151735.1'
    # return 'eva_cannuum_zunla1ref10', 'GCA_000710875.1'
    # return 'eva_lsalmonis_lsalatlcanadafemalev1', 'GCA_000181255.1'


def populate_ids(private_config_xml_file, profile='production', mongo_accession_db='eva_accession_sharded'):
    # get species db name and assembly from evapro
    db_name, assembly = get_db_name_and_assembly_accession(private_config_xml_file)

    # query variants in variant warehouse
    #   check the ids field is empty
    with pymongo.MongoClient(get_mongo_uri_for_eva_profile(profile, private_config_xml_file)) as mongo_handle:
        variants_collection = mongo_handle[db_name]["variants_2_0"]
        projection = {"_id": 1, "files.sid": 1, "chr": 1, "start": 1, "ref": 1, "alt": 1, "type": 1}
        cursor = variants_collection.find({}, projection=projection, limit=2)
        for variant in cursor:
            variant_id, sve_hash = get_ids(assembly, variant)
            hash_variant_id = {sve_hash: variant_id}
            print(hash_variant_id)

            # Get SS ID from accessioning DB
            sve_collection = mongo_handle[mongo_accession_db]['submittedVariantEntity']
            sve_filter = {"_id": sve_hash}
            sve_projection = {"_id": 1, "accession": 1, "rs": 1}
            cursor_accessioning = sve_collection.find_one(sve_filter, projection=sve_projection)

            hash_accession = {sve_hash: cursor_accessioning['accession']}
            # TODO: Support RS IDs
            # if cursor_accessioning.get('rs'):
            #     hash_accession[sve_hash] = f"{hash_accession[sve_hash]}, rs{cursor_accessioning['rs']}"
            print(hash_accession)

            # Insert ids to variant warehouse
            update_statements = []
            for ss_hash, accession in hash_accession.items():
                query = {"_id": hash_variant_id[ss_hash]}
                update = {"$addToSet": {"ids": f"ss{accession}"}}
                update_statements.append(pymongo.UpdateOne(query, update))
            result_insert = variants_collection.bulk_write(requests=update_statements, ordered=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get stats from variant warehouse', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    args = parser.parse_args()
    populate_ids(args.private_config_xml_file)
