import argparse
import base64
import sys
import urllib.request

import psycopg2
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile, get_properties_from_xml_file
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query
from retry import retry

logger = logging_config.get_logger(__name__)


def get_assemblies_from_evapro(private_config_xml_file):
    metadata_handle = psycopg2.connect(get_pg_metadata_uri_for_eva_profile(
        "development", private_config_xml_file), user="evadev")
    query = "select assembly_accession from accessioned_assembly where assembly_accession like 'GCA%'"
    evapro_assemblies = get_all_results_for_query(metadata_handle, query)
    return [asm[0] for asm in evapro_assemblies]


def get_contig_alias_auth(private_config_xml_file):
    properties = get_properties_from_xml_file("production", private_config_xml_file)
    contig_alias_db_username = str(properties['contig-alias.admin-user'])
    contig_alias_db_password = str(properties['contig-alias.admin-password'])

    auth = f"{contig_alias_db_username}:{contig_alias_db_password}"
    encode_auth = auth.encode('ascii')
    base64_auth = base64.b64encode(encode_auth)
    return base64_auth.decode()


@retry(tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def call_admin_endpoint(assemblies, admin_credentials):
    for assembly in assemblies:
        url = f"https://www.ebi.ac.uk/eva/webservices/contig-alias/v1/admin/assemblies/{assembly}"
        headers = {'Authorization': 'Basic ' + admin_credentials}
        request = urllib.request.Request(url, None, headers)
        try:
            with urllib.request.urlopen(request) as response:
                print(f"Assembly {assembly} loaded into contig alias database")
                logger.info(f"Assembly {assembly} loaded into contig alias database")
        except Exception:
            logger.error(f"Error loading assembly {assembly}")


def load_data_to_contig_alias(private_config_xml_file, assembly_list):
    assemblies = assembly_list if assembly_list else get_assemblies_from_evapro(private_config_xml_file)
    admin_credentials = get_contig_alias_auth(private_config_xml_file)
    print(f"Assemblies to be loaded into contig alias database: {assemblies}")
    logger.info(f"Assemblies to be loaded into contig alias database: {assemblies}")
    call_admin_endpoint(assemblies, admin_credentials)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get stats from variant warehouse', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--assembly-list", help="Assembly list e.g. GCA_000181335.4", required=False, nargs='+')
    args = {}
    try:
        args = parser.parse_args()
        load_data_to_contig_alias(args.private_config_xml_file, args.assembly_list)
    except Exception as ex:
        logger.exception(ex)
        sys.exit(1)
    sys.exit(0)
