# Copyright 2020 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This script attempts to get assemblies where a given Contig (genbank accession)
# is present and insert it into a table

from Bio import Entrez
from ebi_eva_common_pyutils.assembly import NCBIAssembly
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.pg_utils import execute_query
from typing import List

import argparse
import os
import psycopg2
import psycopg2.extras
import sys
import wget

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()

Entrez.email = "eva-dev@ebi.ac.uk"
Entrez.tool = "eva"

contig_analysis_table_name = "eva_tasks.eva2469_contig_analysis"
possible_assemblies_table_name = "eva_tasks.eva2469_possible_assemblies_for_contigs"


def create_table_to_collect_possible_assemblies(private_config_xml_file):
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file), user="evadev") \
            as metadata_connection_handle:
        create_table_to_store_possible_assemblies_query = "create table if not exists {0} " \
                                                         "(genbank_accession text, assembly_accession text, " \
                                                          "primary key (genbank_accession, assembly_accession))"\
            .format(possible_assemblies_table_name)

        execute_query(metadata_connection_handle, create_table_to_store_possible_assemblies_query)


def insert_possible_assemblies_for_contig(metadata_connection_handle, contig_accession, possible_assemblies):
    if len(possible_assemblies) > 0:
        with metadata_connection_handle.cursor() as cursor:
            for assembly_accession in possible_assemblies:
                logger.info(f"Inserting possible assembly {assembly_accession} for {contig_accession}...")
                psycopg2.extras.execute_values(cursor,
                                               "INSERT INTO {0} (genbank_accession, assembly_accession) "
                                               "VALUES %s"
                                               .format(possible_assemblies_table_name),
                                               [(contig_accession, assembly_accession)], page_size=100)
            metadata_connection_handle.commit()


def _does_contig_exist_in_assembly(contig_accession: str, assembly_accession: str):
    logger.info(f"Obtaining assembly report for {assembly_accession}...")
    asm = NCBIAssembly(assembly_accession, species_scientific_name=None,
                       reference_directory=None)
    try:
        assembly_report_file_name = os.path.basename(asm.assembly_report_url)
        os.system("rm -f " + assembly_report_file_name)
        wget.download(asm.assembly_report_url)
        output = run_command_with_output(f"Checking if contig {contig_accession} exists in assembly {assembly_accession}",
                                         f'grep -w "{contig_accession}" "{assembly_report_file_name}" | cat',
                                         return_process_output=True)
        return output.strip() != ""
    except Exception as ex:
        logger.error(f"Could not download assembly report for {assembly_accession} due to: " + ex.__str__())
        return False


def _get_assembly_accession_for_id(assembly_id: str):
    assembly_results_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
    assembly_results = Entrez.read(assembly_results_handle)
    return assembly_results["DocumentSummarySet"]["DocumentSummary"][0]["Synonym"]["Genbank"]


def get_assemblies_where_contig_appears(contig_accession: str) -> List[str]:
    contig_results_handle = Entrez.esummary(db="nuccore", id=contig_accession)
    contig_results = list(Entrez.parse(contig_results_handle))
    if len(contig_results) == 0:
        logger.error(f"No records returned for contig {contig_accession} when querying with eutils!")
        return []
    if len(contig_results) > 1:
        logger.error(f"More than one record returned for contig {contig_accession} when querying with eutils!")
        return []

    taxonomy_id = int.__repr__(contig_results[0]['TaxId'])
    assembly_results = Entrez.read(Entrez.esearch(db="assembly", term=f"txid{taxonomy_id}", retmax=100000))
    assemblies_to_check = filter(lambda accession: accession.strip() != "",  [_get_assembly_accession_for_id(assembly_id)
                                  for assembly_id in assembly_results["IdList"]])
    return list(filter(lambda assembly_accession: _does_contig_exist_in_assembly(contig_accession, assembly_accession),
                       assemblies_to_check))


def main():
    parser = argparse.ArgumentParser(description='Get possible assemblies where given Genbank contigs are present',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--private-config-xml-file",
                        help="Full path to private configuration file (ex: /path/to/settings.xml)", required=True)
    parser.add_argument("--eutils-api-key", help="EUtils API key", required=True)
    args = parser.parse_args()

    Entrez.api_key = args.eutils_api_key
    create_table_to_collect_possible_assemblies(args.private_config_xml_file)
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", args.private_config_xml_file),
                          user="evadev") \
            as metadata_connection_handle:
        for contig_accession in sys.stdin:
            contig_accession = contig_accession.strip()
            logger.info(f"Getting possible assemblies for {contig_accession} from EUtils...")
            possible_assemblies = get_assemblies_where_contig_appears(contig_accession)
            insert_possible_assemblies_for_contig(metadata_connection_handle, contig_accession, possible_assemblies)


if __name__ == "__main__":
    main()
