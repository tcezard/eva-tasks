import logging
import re
from argparse import ArgumentParser

import requests
from ebi_eva_common_pyutils.ena_utils import download_xml_from_ena
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.network_utils import json_request
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query

logging.basicConfig(level=logging.INFO)


def get_scientific_name_from_ensembl(taxonomy_id):
    ENSEMBL_REST_API_URL = "https://rest.ensembl.org/taxonomy/id/{0}?content-type=application/json".format(taxonomy_id)
    response = json_request(ENSEMBL_REST_API_URL)
    if 'scientific_name' not in response:
        logging.warning(f'Scientific name not available in Ensembl for taxonomy: {taxonomy_id}')
        return 'Not Available'
    else:
        return response['scientific_name']


def get_scientific_name_from_ena(taxonomy_id):
    xml_root = download_xml_from_ena(f'https://www.ebi.ac.uk/ena/browser/api/xml/{taxonomy_id}')
    xml_taxon = xml_root.xpath('/TAXON_SET/taxon')
    if len(xml_taxon) == 0:
        logging.warning(f'Scientific name not available in ENA for taxonomy: {taxonomy_id}')
        return 'Not Available'
    else:
        return xml_taxon[0].get('scientificName')


def get_scientific_name_from_ncbi(taxonomy_id):
    payload = {'db': 'Taxonomy', 'id': taxonomy_id}
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params=payload)
    match = re.search('<ScientificName>(.+?)</ScientificName>', r.text, re.MULTILINE)
    if match:
        return match.group(1)
    else:
        logging.warning(f'Scientific name not available in ENA for taxonomy: {taxonomy_id}')
        return 'Not Available'


def get_list_of_taxonomy_from_evapro(private_config_xml_file):
    taxonomy_list = []
    query = f"select taxonomy_id from evapro.taxonomy "
    with get_metadata_connection_handle('production_processing', private_config_xml_file) as pg_conn:
        for taxonomy in get_all_results_for_query(pg_conn, query):
            taxonomy_list.append(taxonomy[0])

    return taxonomy_list


def get_scientific_name_from_evapro(private_config_xml_file, taxonomy_id):
    query = f"select scientific_name from evapro.taxonomy " \
            f"where taxonomy_id = {taxonomy_id}"
    with get_metadata_connection_handle('production_processing', private_config_xml_file) as db_conn:
        results = get_all_results_for_query(db_conn, query)
        return results[0][0]


def check_if_name_is_different(taxonomy_list, evapro_names, ensembl_names, ena_names, ncbi_names):
    for taxonomy in taxonomy_list:
        if evapro_names[taxonomy] == ensembl_names[taxonomy] == ena_names[taxonomy] == ncbi_names[taxonomy]:
            pass
        else:
            logging.error(f"""Scientific names differ for taxonomy : {taxonomy}
                EVA : {evapro_names[taxonomy]}
                Ensembl: {ensembl_names[taxonomy]}
                ENA: {ena_names[taxonomy]}
                NCBI: {ncbi_names[taxonomy]}
            """)


def main():
    argparser = ArgumentParser(description='Retrieve scientific names form different sources like evparo, Ensembl, '
                                           'ENA and NCBI and check if there is any difference between them')
    argparser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml",
                           required=True)
    args = argparser.parse_args()

    taxonomy_list = get_list_of_taxonomy_from_evapro(args.private_config_xml_file)

    evapro_names = {}
    ensembl_names = {}
    ena_names = {}
    ncbi_names = {}
    count = 1
    for taxonomy in taxonomy_list:
        logging.info(f'{count}. {taxonomy}')

        evapro_names[taxonomy] = get_scientific_name_from_evapro(args.private_config_xml_file, taxonomy)
        ensembl_names[taxonomy] = get_scientific_name_from_ensembl(taxonomy)
        ena_names[taxonomy] = get_scientific_name_from_ena(taxonomy)
        ncbi_names[taxonomy] = get_scientific_name_from_ncbi(taxonomy)

        count = count + 1

    check_if_name_is_different(taxonomy_list, evapro_names, ensembl_names, ena_names, ncbi_names)


if __name__ == "__main__":
    main()
