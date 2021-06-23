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

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.assembly import NCBIAssembly
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.pg_utils import execute_query

import argparse
import os
import psycopg2
import psycopg2.extras
import sys
import traceback
import wget

# This script is adapted from  https://github.com/EBIvariation/eva-tasks/blob/5411bd9a1187140480efd04fa159bd0b4aad0123/tasks/eva_2150/collect_assembly_report_genbank_contigs.py
# Few changes have been made to add 2 additional columns to the table:
# 1) RefSeq accession and
# 2) RefSeq-Genbank contig equivalence relationship
logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()
asm_report_contigs_table_name = "eva_tasks.eva2469_asm_report_genbank_contigs"


def create_table_to_collect_assembly_report_genbank_contigs(private_config_xml_file):    
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file), user="evadev") \
            as metadata_connection_handle:
        create_table_to_store_asm_report_contigs_query = "create table if not exists {0} " \
                                                         "(assembly_accession text, genbank_accession text, " \
                                                         "chromosome_name text, " \
                                                         "is_equivalent_genbank_available boolean, " \
                                                         "refseq_accession text)".format(asm_report_contigs_table_name)

        execute_query(metadata_connection_handle, create_table_to_store_asm_report_contigs_query)
        

def insert_contigs_to_db(metadata_connection_handle, contig_info_list):
    if len(contig_info_list) > 0:
        with metadata_connection_handle.cursor() as cursor:
            psycopg2.extras.execute_values(cursor,
                                           "INSERT INTO {0} (assembly_accession, genbank_accession, chromosome_name, "
                                           "is_equivalent_genbank_available, refseq_accession) "
                                           "VALUES %s".format(asm_report_contigs_table_name), contig_info_list,
                                           page_size=100)


def collect_assembly_report_genbank_contigs(private_config_xml_file, assembly_accession):
    try:
        with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                              user="evadev") \
                as metadata_connection_handle:
            asm = NCBIAssembly(assembly_accession, species_scientific_name=None,
                               reference_directory=None)
            assembly_report_file_name = os.path.basename(asm.assembly_report_url)
            os.system("rm -f " + assembly_report_file_name)
            wget.download(asm.assembly_report_url)

            insert_chunk_size = 100
            contig_info_list = []
            for line in open(assembly_report_file_name, 'r'):
                if not line.strip().startswith("#"):
                    line_components = line.strip().split("\t")
                    chromosome_name, genbank_accession, accession_equivalence, refseq_accession = \
                        line_components[0], line_components[4], line_components[5], line_components[6]
                    # Equivalence "Relationship" column in the assembly report indicates if
                    # Genbank and RefSeq contig accessions are equivalent
                    is_equivalent_genbank_available = (accession_equivalence.strip() == "=")
                    contig_info_list.append((assembly_accession, genbank_accession, chromosome_name,
                                             is_equivalent_genbank_available, refseq_accession))
                    if len(contig_info_list) == insert_chunk_size:
                        insert_contigs_to_db(metadata_connection_handle, contig_info_list)
                        contig_info_list = []
            insert_contigs_to_db(metadata_connection_handle, contig_info_list)
    except Exception:
        logger.error(traceback.format_exc())


def main():
    parser = argparse.ArgumentParser(description='Collect assembly report Genbank contigs for a given assembly',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--private-config-xml-file",
                        help="Full path to private configuration file (ex: /path/to/settings.xml)", required=True)
    args = parser.parse_args()

    create_table_to_collect_assembly_report_genbank_contigs(args.private_config_xml_file)
    for assembly_accession in sys.stdin:
        collect_assembly_report_genbank_contigs(args.private_config_xml_file, assembly_accession.strip())


if __name__ == "__main__":
    main()
