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

from eva_2150 import init_logger
from ebi_eva_common_pyutils.variation import assembly_utils
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.pg_utils import execute_query

import click
import os
import psycopg2
import psycopg2.extras
import sys
import traceback
import wget


logger = init_logger()
asm_report_contigs_table_name = "eva_tasks.eva2150_asm_report_genbank_contigs"


def create_table_to_collect_assembly_report_genbank_contigs(private_config_xml_file):    
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file), user="evadev") \
            as metadata_connection_handle:
        create_table_to_store_asm_report_contigs_query = "create table if not exists {0} " \
                                                         "(assembly_accession text, contig_accession text, " \
                                                         "chromosome_name text)".format(asm_report_contigs_table_name)

        execute_query(metadata_connection_handle, create_table_to_store_asm_report_contigs_query)
        

def insert_contigs_to_db(metadata_connection_handle, contig_info_list):
    if len(contig_info_list) > 0:
        with metadata_connection_handle.cursor() as cursor:
            psycopg2.extras.execute_values(cursor,
                                           "INSERT INTO {0} (assembly_accession, contig_accession, chromosome_name) "
                                           "VALUES %s".format(asm_report_contigs_table_name), contig_info_list,
                                           page_size=100)


def collect_assembly_report_genbank_contigs(private_config_xml_file, assembly_accession):
    try:
        with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                              user="evadev") \
                as metadata_connection_handle:
            assembly_report_url = assembly_utils.get_assembly_report_url(assembly_accession)
            assembly_report_file_name = os.path.basename(assembly_report_url)
            os.system("rm -f " + assembly_report_file_name)
            wget.download(assembly_report_url)

            insert_chunk_size = 100
            contig_info_list = []
            for line in open(assembly_report_file_name, 'r'):
                if not line.strip().startswith("#"):
                    line_components = line.strip().split("\t")
                    chromosome_name, genbank_accession = line_components[0], line_components[4]
                    contig_info_list.append((assembly_accession, genbank_accession, chromosome_name))
                    if len(contig_info_list) == insert_chunk_size:
                        insert_contigs_to_db(metadata_connection_handle, contig_info_list)
                        contig_info_list = []
            insert_contigs_to_db(metadata_connection_handle, contig_info_list)
    except Exception:
        logger.error(traceback.format_exc())


@click.option("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
@click.command()
def main(private_config_xml_file):
    create_table_to_collect_assembly_report_genbank_contigs(private_config_xml_file)
    for assembly_accession in sys.stdin:
        collect_assembly_report_genbank_contigs(private_config_xml_file, assembly_accession.strip())


if __name__ == "__main__":
    main()
