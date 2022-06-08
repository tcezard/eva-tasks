import argparse
import base64
import urllib.request

from ebi_eva_common_pyutils.config_utils import get_properties_from_xml_file
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query
from retry import retry

logger = logging_config.get_logger(__name__)


def get_assemblies_from_evapro(private_config_xml_file):
    with get_metadata_connection_handle("development", private_config_xml_file) as pg_conn:
        query = "select assembly_accession from accessioned_assembly where assembly_accession like 'GCA%'"
        evapro_assemblies = get_all_results_for_query(pg_conn, query)
        return [asm[0] for asm in evapro_assemblies]


def get_contig_alias_auth(private_config_xml_file):
    properties = get_properties_from_xml_file("development", private_config_xml_file)
    contig_alias_username = str(properties['contig-alias.admin-user'])
    contig_alias_password = str(properties['contig-alias.admin-password'])

    auth = f"{contig_alias_username}:{contig_alias_password}"
    encode_auth = auth.encode('ascii')
    base64_auth = base64.b64encode(encode_auth)
    return base64_auth.decode()


def load_assembly_to_contig_alias(assemblies, admin_credentials, delete_insert):
    for assembly in assemblies:
        url = f"https://wwwdev.ebi.ac.uk/eva/webservices/contig-alias/v1/admin/assemblies/{assembly}"
        headers = {'Authorization': 'Basic ' + admin_credentials}
        if delete_insert == 'True':
            #delete request for assembly
            del_request(assembly, url, headers)

        # insert request for assembly
        insert_request(assembly, url, headers)


def del_request(assembly, url, headers):
    try:
        create_and_execute_request(url, 'DELETE', headers)
        logger.warning(f'Assembly info Deleted successfully for Assembly {assembly}')
    except Exception as e:
        logger.error(f'Error Deleting assembly {assembly}')


def insert_request(assembly, url, headers):
    try:
        create_and_execute_request(url, 'PUT', headers)
        logger.warning(f'Assembly info inserted successfully for Assembly {assembly}')
    except Exception as e:
        logger.error(f'Error inserting assembly {assembly}')


#@retry(tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def create_and_execute_request(url, method, headers):
    request = urllib.request.Request(url, None, headers, method=method)
    urllib.request.urlopen(request)


def load_data_to_contig_alias(private_config_xml_file, assembly_list, delete_insert):
    assemblies = assembly_list if assembly_list else get_assemblies_from_evapro(private_config_xml_file)
    admin_credentials = get_contig_alias_auth(private_config_xml_file)
    logger.info(f"Assemblies to be loaded into contig alias database: {assemblies}")
    load_assembly_to_contig_alias(assemblies, admin_credentials, delete_insert)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Load data into contig alias database', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--assembly-list", help="Assembly list e.g. GCA_000181335.4", required=False, nargs='+')
    parser.add_argument("--delete-insert", choices=('True', 'False'), help="Delete and re insert Assembly information",
                        required=False, default='False')
    args = parser.parse_args()

    load_data_to_contig_alias(args.private_config_xml_file, args.assembly_list, args.delete_insert)
