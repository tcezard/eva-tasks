import argparse
import os

import requests
from ebi_eva_common_pyutils.config_utils import get_contig_alias_db_creds_for_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query
from retry import retry

logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)


class InternalServerError(Exception):
    pass


def get_assemblies_from_evapro(profile, private_config_xml_file):
    with get_metadata_connection_handle(profile, private_config_xml_file) as pg_conn:
        query = "select distinct assembly_accession from evapro.accessioned_assembly where assembly_accession like 'GCA%'" \
                " union " \
                "select distinct assembly_accession from eva_progress_tracker.remapping_tracker where assembly_accession like 'GCA%'"
        evapro_assemblies = get_all_results_for_query(pg_conn, query)
        return [asm[0] for asm in evapro_assemblies]


def load_assembly_to_contig_alias(assemblies, contig_alias_url, contig_alias_user, contig_alias_pass, overwrite):
    logger.info(f"A total of {len(assemblies)} assemblies to be loaded into contig-alias database: {assemblies}")

    for assembly in assemblies:
        full_url = os.path.join(contig_alias_url, f'v1/admin/assemblies/{assembly}')
        if overwrite:
            # delete request for assembly
            del_request(assembly, full_url, contig_alias_user, contig_alias_pass)

        # insert request for assembly
        insert_request(assembly, full_url, contig_alias_user, contig_alias_pass)


@retry(InternalServerError, tries=3, delay=2, backoff=1.5, jitter=(1, 3))
def del_request(assembly, url, user, password):
    response = requests.delete(url, auth=(user, password))
    if response.status_code == 200:
        logger.info(f'Assembly accession {assembly} successfully deleted from Contig-Alias DB')
    elif response.status_code == 500:
        logger.error(f'Assembly accession {assembly} could not be deleted. Response: {response.text}')
        raise InternalServerError
    else:
        logger.error(f'Assembly accession {assembly} could not be deleted. Response: {response.text}')


@retry(InternalServerError, tries=10, delay=2, backoff=1.5, jitter=(1, 3))
def insert_request(assembly, url, user, password):
    response = requests.put(url, auth=(user, password))
    if response.status_code == 200:
        logger.info(f'Assembly accession {assembly} successfully added to Contig-Alias DB')
    elif response.status_code == 409:
        logger.warning(f'Assembly accession {assembly} already exist in Contig-Alias DB. Response: {response.text}')
    elif response.status_code == 500:
        logger.error(f'Could not save Assembly accession {assembly} to Contig-Alias DB. Error : {response.text}')
        raise InternalServerError
    else:
        logger.error(f'Could not save Assembly accession {assembly} to Contig-Alias DB. Error : {response.text}')


def load_data_to_contig_alias(private_config_xml_file, profile, assembly_list, overwrite):
    assemblies = assembly_list if assembly_list else get_assemblies_from_evapro(profile, private_config_xml_file)
    contig_alias_url, contig_alias_user, contig_alias_pass = get_contig_alias_db_creds_for_profile(
        profile, private_config_xml_file)

    load_assembly_to_contig_alias(assemblies, contig_alias_url, contig_alias_user, contig_alias_pass, overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Load data into contig alias database', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--profile", choices=('localhost', 'development', 'production'),
                        help="Profile to decide whether to run this for development or production", required=True)
    parser.add_argument("--assembly-list", help="Assembly list e.g. GCA_000181335.4", required=False, nargs='+')
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Whether to delete and re-insert assembly information")

    args = parser.parse_args()

    load_data_to_contig_alias(args.private_config_xml_file, args.profile, args.assembly_list, args.overwrite)

