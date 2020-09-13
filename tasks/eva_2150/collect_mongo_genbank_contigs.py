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
from ebi_eva_common_pyutils.variation import contig_utils
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile, get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.pg_utils import execute_query, get_all_results_for_query

import click
import psycopg2
import psycopg2.extras
import sys
import traceback

from pymongo import MongoClient

logger = init_logger()
mongo_genbank_contigs_table_name = "eva_tasks.eva2150_mongo_genbank_contigs"


def get_chromosome_name_from_asm_report(metadata_connection_handle, assembly_accession, contig_accession):
    query = "select chromosome_name from eva_tasks.eva2150_asm_report_genbank_contigs " \
            "where assembly_accession = '{0}' and contig_accession = '{1}'"\
        .format(assembly_accession, contig_accession)
    results = get_all_results_for_query(metadata_connection_handle, query)
    if len(results) == 0:
        return None
    if len(results) > 1:
        logger.error(
            "More than one chromosome name found for assembly: {0} and contig: {1}".format(assembly_accession,
                                                                                           contig_accession))
    return results[0][0]


def create_table_to_collect_mongo_genbank_contigs(private_config_xml_file):
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file), user="evadev") \
            as metadata_connection_handle:
        create_table_to_store_asm_report_contigs_query = "create table if not exists {0} " \
                                                         "(source text, assembly_accession text, " \
                                                         "study text, contig_accession text, chromosome_name text, " \
                                                         "num_entries_in_db bigint, is_contig_in_asm_report boolean, " \
                                                         "primary key(source, assembly_accession, study, " \
                                                         "contig_accession))"\
            .format(mongo_genbank_contigs_table_name)

        execute_query(metadata_connection_handle, create_table_to_store_asm_report_contigs_query)


def insert_contigs_to_db(metadata_connection_handle, contig_info_list):
    if len(contig_info_list) > 0:
        with metadata_connection_handle.cursor() as cursor:
            psycopg2.extras.execute_values(cursor,
                                           "INSERT INTO {0} "
                                           "(source, assembly_accession, study, contig_accession, chromosome_name, "
                                           "num_entries_in_db, is_contig_in_asm_report) "
                                           "VALUES %s".format(mongo_genbank_contigs_table_name), contig_info_list,
                                           page_size=100)


def insert_contig_info_to_db(collection, assembly_accession, metadata_connection_handle, mongo_connection_handle,
                             assembly_attribute_prefix=""):
    collection_handle = mongo_connection_handle["eva_accession_sharded"][collection]
    with collection_handle.aggregate([{'$match': {assembly_attribute_prefix + 'seq': assembly_accession}},
                                      {'$group': {'_id': {'study': '$' + assembly_attribute_prefix + 'study',
                                                          'contig': '$' + assembly_attribute_prefix + 'contig'},
                                                  'count': {'$sum': 1}}},
                                      {"$project": {"study": "$_id.study", "contig": "$_id.contig",
                                                    "count": 1, "_id": 0}}
                                      ], allowDiskUse=True) as cursor:
        insert_chunk_size = 100
        contig_info_list = []
        for result in cursor:
            is_contig_in_asm_report = False
            study = result["study"][0] if assembly_attribute_prefix else result["study"]
            genbank_accession = result["contig"][0] if assembly_attribute_prefix else result["contig"]
            count = result["count"]
            chromosome_name = get_chromosome_name_from_asm_report(metadata_connection_handle,
                                                                  assembly_accession, genbank_accession)
            if chromosome_name is not None:
                is_contig_in_asm_report = True
            else:
                chromosome_name = contig_utils.get_chromosome_name_for_contig_accession(genbank_accession)
            contig_info_list.append((collection, assembly_accession, study, genbank_accession,
                                     chromosome_name, count, is_contig_in_asm_report))
            if len(contig_info_list) == insert_chunk_size:
                insert_contigs_to_db(metadata_connection_handle, contig_info_list)
                contig_info_list = []
        insert_contigs_to_db(metadata_connection_handle, contig_info_list)


def collect_mongo_genbank_contigs(private_config_xml_file, assembly_accession):
    try:
        with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                              user="evadev") \
                as metadata_connection_handle, MongoClient(get_mongo_uri_for_eva_profile("development",
                                                                                         private_config_xml_file)) \
                as mongo_connection_handle:
            main_collections = ["dbsnpSubmittedVariantEntity", "submittedVariantEntity"]
            for collection in main_collections:
                insert_contig_info_to_db(collection, assembly_accession,
                                         metadata_connection_handle, mongo_connection_handle)
            ops_collections = ["dbsnpSubmittedVariantOperationEntity", "submittedVariantOperationEntity"]
            for collection in ops_collections:
                insert_contig_info_to_db(collection, assembly_accession,
                                         metadata_connection_handle, mongo_connection_handle,
                                         assembly_attribute_prefix="inactiveObjects.")
    except Exception:
        logger.error(traceback.format_exc())


@click.option("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
@click.command()
def main(private_config_xml_file):
    create_table_to_collect_mongo_genbank_contigs(private_config_xml_file)
    for assembly_accession in sys.stdin:
        logger.info("Processing assembly: " + assembly_accession)
        collect_mongo_genbank_contigs(private_config_xml_file, assembly_accession.strip())


if __name__ == "__main__":
    main()
