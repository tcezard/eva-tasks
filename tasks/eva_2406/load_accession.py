#!/usr/bin/env python

from argparse import ArgumentParser


import psycopg2
import psycopg2.extras
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.pg_utils import execute_query, get_pg_connection_handle, get_all_results_for_query

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()

"""Simple script to load the summarise output of submitted variants collection in accessioning database."""


def parse_accessioning_counts_from_reports(accession_report_counts_file):
    accession_count = {}
    with open(accession_report_counts_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split('\t')
            accession_count[sp_line[0]] = int(sp_line[1])
    return accession_count


def parse_accessioning_counts_from_db(accession_db_counts_file):
    """Parse file that contains the accessioning counts.
    Expects 4 columns: assembly accession, taxonomy id, project id, datestamp, nb_ss_variants """
    accession_count = []
    with open(accession_db_counts_file) as open_file:
        for line in open_file:
            accession_count.append(line.strip().split('\t'))
    return accession_count


def create_table_accession_counts(private_config_xml_file):
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                          user="evadev") as metadata_connection_handle:
        query_create_table = (
            'CREATE TABLE IF NOT EXISTS eva_stats.submitted_variants_load_counts '
            '(source TEXT, taxid INTEGER, assembly_accession TEXT, project_accession TEXT, date_loaded TIMESTAMP, '
            'number_submitted_variants BIGINT NOT NULL, '
            'primary key(taxid, assembly_accession, project_accession, date_loaded))'
        )
    execute_query(metadata_connection_handle, query_create_table)


def insert_accession_counts_to_db(private_config_xml_file, accession_counts, source):
    if len(accession_counts) > 0:
        with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file),
                              user="evadev") as metadata_connection_handle:
            with metadata_connection_handle.cursor() as cursor:
                query_insert = (
                    'INSERT INTO eva_stats.submitted_variants_load_counts '
                    '(source, assembly_accession, taxid, project_accession, date_loaded, number_submitted_variants) '
                    'VALUES %s '
                    'ON CONFLICT (taxid, assembly_accession, project_accession, date_loaded) '
                    'DO UPDATE SET number_submitted_variants = EXCLUDED.number_submitted_variants'
                )
                psycopg2.extras.execute_values(cursor, query_insert, accession_counts, ("('" + source + "', %s, %s, %s, %s, %s)"))


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--private_config_xml_file', help='Path to the file containing the ', required=True)
    argparse.add_argument('--accession_counts', help='Path to the file that contain counts per project id', required=True)
    argparse.add_argument('--source', help='Source of the study', choices=['EVA', 'dbSNP'], required=True)
    args = argparse.parse_args()

    accession_count = parse_accessioning_counts_from_db(args.accession_counts)
    create_table_accession_counts(args.private_config_xml_file)
    insert_accession_counts_to_db(args.private_config_xml_file, accession_count, args.source)


if __name__ == "__main__":
    main()
