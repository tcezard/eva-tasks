import hashlib
import argparse
import pymongo
import psycopg2
import requests
import traceback
import logging
from csv import DictReader, excel_tab
from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile


def generate_update_statement(variant_id, ss_accession, rs_accession):
    query = {"_id": variant_id}
    if rs_accession:
        update = {"$addToSet": {"ids": {"$each": [f"ss{ss_accession}", f"rs{rs_accession}"]}}}
    else:
        update = {"$addToSet": {"ids": f"ss{ss_accession}"}}
    return pymongo.UpdateOne(query, update)


def get_from_accessioning_db(mongo_handle, mongo_accession_db, sve_hash):
    # Get SS ID from accessioning DB
    sve_collection = mongo_handle[mongo_accession_db]['submittedVariantEntity']
    sve_filter = {"_id": sve_hash}
    sve_projection = {"_id": 1, "accession": 1, "rs": 1}
    cursor_accessioning = sve_collection.find_one(sve_filter, projection=sve_projection)
    ss_accession = cursor_accessioning['accession']
    # Get RS ID is variant is clustered
    rs_accession = cursor_accessioning.get('rs')
    return ss_accession, rs_accession


def get_variant_warehouse_cursor(variants_collection):
    projection = {"_id": 1, "files.sid": 1, "chr": 1, "start": 1, "ref": 1, "alt": 1, "type": 1}
    return variants_collection.find({}, projection=projection)


def get_assembly_report_rows(assembly_report_path):
    with open(assembly_report_path) as open_file:
        headers = None
        # Parse the assembly report file to find the header then stop
        for line in open_file:
            if line.lower().startswith("# sequence-name") and "sequence-role" in line.lower():
                headers = line.strip().split('\t')
                break
        reader = DictReader(open_file, fieldnames=headers, dialect=excel_tab)
        for record in reader:
            yield record


def get_genbank_contig_alias(contig):
    url = 'http://ves-hx-de.ebi.ac.uk:8080/eva/webservices/contig-alias/v1/chromosomes/refseq/' + contig
    response = requests.get(url)
    genbank = response.json()['_embedded']['chromosomeEntities'][0]['genbank']
    return genbank


def get_SHA1(variant_id):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    return hashlib.sha1(variant_id.encode()).hexdigest().upper()


def get_ids(assembly, variant_query_result):
    """Get the id fields as a pair (id_variant_warehouse, sve_hash)"""
    ss_hash_variant_id = {}
    for file in variant_query_result['files']:
        study = file['sid']
        genbank_chr = get_genbank_contig_alias(variant_query_result['chr'])
        start = variant_query_result['start']
        ref = variant_query_result['ref']
        alt = variant_query_result['alt']
        id_variant_warehouse = variant_query_result['_id']
        sve_hash = get_SHA1(f"{assembly}_{study}_{genbank_chr}_{start}_{ref}_{alt}")
        ss_hash_variant_id[sve_hash] = id_variant_warehouse
    return ss_hash_variant_id


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
    with pymongo.MongoClient(get_mongo_uri_for_eva_profile(profile, private_config_xml_file)) as mongo_handle:
        variants_collection = mongo_handle[db_name]["variants_2_0"]
        cursor = get_variant_warehouse_cursor(variants_collection)
        update_statements = []
        for variant in cursor:
            for sve_hash, variant_id in get_ids(assembly, variant).items():
                ss_accession, rs_accession = get_from_accessioning_db(mongo_handle, mongo_accession_db, sve_hash)
                update_statements.append(generate_update_statement(variant_id, ss_accession, rs_accession))
        result_insert = variants_collection.bulk_write(requests=update_statements, ordered=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get stats from variant warehouse', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    args = parser.parse_args()
    populate_ids(args.private_config_xml_file)
