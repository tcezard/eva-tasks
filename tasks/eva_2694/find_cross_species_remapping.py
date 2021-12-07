from argparse import ArgumentParser

import psycopg2
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.ncbi_utils import get_ncbi_assembly_dicts_from_term
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query


def main(private_config_xml_file):
    private_config_xml_file
    with psycopg2.connect(get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file), user="evadev") as pg_conn:
        query = (
            "select distinct origin_assembly_accession, assembly_accession "
            "from eva_progress_tracker.remapping_tracker "
            "where origin_assembly_accession!=assembly_accession and num_ss_ids>0"
        )
        source_assembly, target_assembly = get_all_results_for_query(pg_conn, query)
        source_assembly_info = get_ncbi_assembly_dicts_from_term(source_assembly)[0]
        target_assembly_info = get_ncbi_assembly_dicts_from_term(target_assembly)[0]
        source_taxid = source_assembly_info['speciestaxid']
        target_taxid = target_assembly_info['speciestaxid']
        source_organism = source_assembly_info['organism']
        target_organism = target_assembly_info['organism']
        if source_taxid != target_taxid:
            print(f'{source_assembly} and {target_assembly} have different source species {source_organism} != {target_organism}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    args = parser.parse_args()
    main(args.private_config_xml_file)

