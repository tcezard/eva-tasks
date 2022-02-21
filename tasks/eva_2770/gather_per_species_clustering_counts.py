import argparse
import os

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query, execute_query

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


counts_table_name = 'dbsnp_ensembl_species.release_rs_statistics_per_species'
tracker_table_name = 'eva_progress_tracker.clustering_release_tracker'
shell_script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bash')

# TODO other columns are "new" which is just the diff with release 2?
id_to_column = {
    'current': 'current_rs',
    'multimap': 'multi_mapped_rs',
    'merged': 'merged_rs',
    'deprecated': 'deprecated_rs',
    'merged_deprecated': 'merged_deprecated_rs',
    'unmapped': 'unmapped_rs'
}


def write_counts_to_table(private_config_xml_file, counts):
    with get_metadata_connection_handle('development', private_config_xml_file) as db_conn:
        all_columns = ['taxonomy_id', 'scientific_name', 'release_folder', 'release_version']
        all_columns.extend(id_to_column.values())
        all_values = [f"({','.join(species_counts[c] for c in all_columns)})" for species_counts in counts]
        insert_query = f"insert into {counts_table_name} " \
                       f"({','.join(all_columns)}) " \
                       f"values {','.join(all_values)}"
        execute_query(db_conn, insert_query)


def get_scientific_name_and_taxonomy(private_config_xml_file, release_version, species_folder):
    query = f"select taxonomy, scientific_name from {tracker_table_name} " \
            f"where release_version={release_version} " \
            f"and release_folder_name='{species_folder}' " \
            f"and should_be_released"
    with get_metadata_connection_handle('development', private_config_xml_file) as db_conn:
        results = get_all_results_for_query(db_conn, query)
    if len(results) != 1:
        raise ValueError(f'Failed to get scientific name and taxonomy for {species_folder}')
    return results[0][0], results[0][1]


def run_count_script(script_name, species_dir, metric_id):
    run_command_with_output(
        f'Run {script_name}',
        f'bsub {os.path.join(shell_script_dir, script_name)} {species_dir} {metric_id}'
    )
    return f'{os.path.basename(species_dir)}_count_{metric_id}_rsid.log'


def gather_counts(private_config_xml_file, release_version, release_dir):
    results = {}
    for species_dir in os.listdir(release_dir):
        full_species_dir = os.path.abspath(species_dir)
        taxid, sci_name = get_scientific_name_and_taxonomy(private_config_xml_file, release_version, species_dir)
        per_species_results = {
            'taxonomy_id': taxid,
            'scientific_name': sci_name,
            'release_folder': species_dir,
            'release_version': release_version
        }

        for metric_id in id_to_column.keys():
            if metric_id in {'current', 'merged', 'multimap'}:
                output_log = run_count_script('count_rs_for_release.sh', full_species_dir, metric_id)
            elif metric_id in {'deprecated', 'merged_deprecated'}:
                output_log = run_count_script('count_rs_for_release_for_txt.sh', full_species_dir, metric_id)
            else:
                output_log = run_count_script('count_rs_for_release_unmapped.sh', full_species_dir, '')

            with open(output_log) as f:
                total = sum(int(l.strip()) for l in f)
            os.remove(output_log)
            per_species_results[id_to_column[metric_id]] = total

        results[taxid] = per_species_results
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Parse all the release output to get RS statistics per species')
    parser.add_argument("--release_root_path", type=str,
                        help="base directory where all the release was run.", required=True)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--release_version", type=int, help="release version to store in the table", required=True)

    args = parser.parse_args()
    counts = gather_counts(args.private_config_xml_file, args.release_version, args.release_root_path)
    write_counts_to_table(args.private_config_xml_file, counts)


if __name__ == '__main__':
    main()
