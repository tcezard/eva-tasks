import sys
from argparse import ArgumentParser
from csv import DictReader, excel_tab

import psycopg2
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query

logger = log_cfg.get_logger(__name__)


def get_contig_genbank(assembly_report_path):
    genbank_to_row = {}
    with open(assembly_report_path) as open_file:
        headers = None
        # Parse the assembly report file to find the header then stop
        for line in open_file:
            if line.lower().startswith("# sequence-name") and "sequence-role" in line.lower():
                headers = line.strip().split('\t')
                break
        reader = DictReader(open_file, fieldnames=headers, dialect=excel_tab)
        for record in reader:
            genbank_to_row[record['GenBank-Accn']] = record
    return genbank_to_row


def get_contigs_accessions_for(pg_conn, accession):
    db_query = ("select distinct contig_accession "
                "from eva_tasks.eva2150_mongo_genbank_contigs "
                "where source='%s' AND assembly_accession='%s'")

    all_eva_contigs = [contig for contig, in get_all_results_for_query(pg_conn, db_query % ('submittedVariantEntity', accession))]
    all_dbsnp_contigs = [contig for contig, in get_all_results_for_query(pg_conn, db_query % ('dbsnpSubmittedVariantOperationEntity', accession))]
    return all_eva_contigs, all_dbsnp_contigs


def main():
    argparser = ArgumentParser()
    argparser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    argparser.add_argument("--assembly_accession", help="GCA_000003205.1", required=True)
    argparser.add_argument("--assembly_report_path", help="path to the report to check contigs against", required=True)
    args = argparser.parse_args()

    genbank_to_row = get_contig_genbank(args.assembly_report_path)

    log_cfg.add_stdout_handler()

    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", args.private_config_xml_file),  user="evadev") as pg_conn:
        eva_contigs, dbSNP_contigs = get_contigs_accessions_for(pg_conn, args.assembly_accession)

        for contig in eva_contigs:
            if contig not in genbank_to_row:
                logger.warning('For assembly {} contig {} found in EVA is not genbank in the report {}'.format(args.assembly_accession, contig, args.assembly_report_path))
        for contig in dbSNP_contigs:
            if contig not in genbank_to_row:
                logger.warning('For assembly {} contig {} found in dbSNP is not genbank in the report {}'.format(args.assembly_accession, contig, args.assembly_report_path))

    return 0


if __name__ == "__main__":
    sys.exit(main())
