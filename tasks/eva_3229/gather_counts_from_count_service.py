from argparse import ArgumentParser
from collections import defaultdict, Counter
from urllib.parse import urlsplit

import psycopg2
from ebi_eva_common_pyutils.config_utils import get_accession_pg_creds_for_profile
from ebi_eva_common_pyutils.logger import AppLogger

from ebi_eva_common_pyutils.common_utils import pretty_print
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query


class CountStats(AppLogger):

    def __init__(self, profile, settings_file):
        self.profile = profile
        self.settings_file = settings_file

    def get_job_tracker_connection_handle(self):
        url, username, password = get_accession_pg_creds_for_profile(self.profile, self.settings_file)
        return psycopg2.connect(urlsplit(url).path, user=username, password=password)

    def get_clustering_counts_between_dates(self, start_date, end_date):
        # Dates are strings of the form YYYY-MM-DD (inclusive of start, exclusive of end)
        with get_metadata_connection_handle(self.profile, self.settings_file) as pg_conn:
            query = (
                f"SELECT identifier->>'assembly' as assembly, metric, count, timestamp FROM evapro.count_stats "
                f"WHERE process = 'clustering' "
                f"AND timestamp >= '{start_date}' AND timestamp < '{end_date}' "
            )
            results = get_all_results_for_query(pg_conn, query)

        results_by_timestamp = defaultdict(dict)
        for row in results:
            results_by_timestamp[(row[0], row[3])][row[1]] = int(row[2])
        results_by_assembly = defaultdict(Counter)
        for k, v in results_by_timestamp.items():
            assembly = k[0]
            for metric, count in v.items():
                # clustered_variants_created is tricky and needs to be disambiguated
                if metric == 'clustered_variants_created' and count > 0:
                    # If submitted variants are newly clustered in this step, probably a true creation of RS
                    if v['submitted_variants_clustered'] > 0:
                        results_by_assembly[assembly][metric] += count
                    # If clustered variants are split, we've also created a new RS (second clause is due to the fact
                    # that currently the merge step counts creation of split *candidates* under this metric as well)
                    elif v['clustered_variants_rs_split'] > 0 and v['clustered_variants_merge_operations'] == 0:
                        results_by_assembly[assembly][metric] += count
                    # If clustered variants are created and deprecated in a single step, probably a remediation
                    # These can be ignored
                    elif v['clustered_variants_deprecated'] > 0:
                        continue
                    # Otherwise this is probably creation of clustered variants due to remapping
                    else:
                        results_by_assembly[assembly]['clustered_variants_remapped'] += count
                # Only counted when merge step removes clustered variants
                elif metric == 'clustered_variants_updated' and count > 0:
                    results_by_assembly[assembly]['clustered_variants_remapped'] -= count
                else:
                    results_by_assembly[assembly][metric] += count
        self.report(results_by_assembly)

    def get_clustering_counts_between_dates_with_job_tracker(self, start_date, end_date):
        created_rs_steps = ['CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP', 'STUDY_CLUSTERING_STEP', 'PROCESS_RS_SPLIT_CANDIDATES_STEP']
        remapped_rs_steps = ['CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP', 'PROCESS_RS_MERGE_CANDIDATES_STEP']
        deprecated_rs_steps = ['DEPRECATE_STUDY_SUBMITTED_VARIANTS_STEP', 'stepsForEVA3101Deprecation', 'stepsForSSDeprecation']
        all_steps_joined = ",".join(f"'{step}'" for step in created_rs_steps + remapped_rs_steps + deprecated_rs_steps)
        with self.get_job_tracker_connection_handle() as jt_conn:
            jt_query = (
                "SELECT step_name, start_time, end_time, status FROM batch_step_execution "
                f"WHERE start_time >= '{start_date}' AND end_time < '{end_date}' "
                f"AND step_name IN ({all_steps_joined}) "
            )
            jt_results = get_all_results_for_query(jt_conn, jt_query)

        results_by_assembly = defaultdict(Counter)
        # Construct query that keeps track of step name and status for each time range
        step_name_case = ""
        status_case = ""
        for step_name, start_time, end_time, status in jt_results:
            step_name_case += f" WHEN timestamp >= '{start_time}' AND timestamp <= '{end_time}' THEN '{step_name}'"
            status_case += f" WHEN timestamp >= '{start_time}' AND timestamp <= '{end_time}' THEN '{status}'"
        query = (
            f"SELECT identifier->>'assembly' as assembly, metric, count, "
            f"CASE {step_name_case} END step_name, "
            f"CASE {status_case} END status "
            f"FROM evapro.count_stats WHERE process = 'clustering'"
        )
        with get_metadata_connection_handle(self.profile, self.settings_file) as pg_conn:
            results = get_all_results_for_query(pg_conn, query)
            for asm, metric, count, step_name, status in results:
                # Follows same logic as get_clustering_counts_between_dates to disambiguate metrics, just this time
                # using the known step names
                if metric == 'clustered_variants_created' and count > 0:
                    if step_name in created_rs_steps:
                        results_by_assembly[asm][metric] += count
                    elif step_name in remapped_rs_steps:
                        results_by_assembly[asm]['clustered_variants_remapped'] += count
                    else:
                        continue
                elif metric == 'clustered_variants_updated' and count > 0:
                    results_by_assembly[asm]['clustered_variants_remapped'] -= count
                else:
                    results_by_assembly[asm][metric] += count

        self.report(results_by_assembly)
        # Dump raw results for debugging purposes
        with open('results.tsv', 'w+') as f:
            f.write('\n'.join('\t'.join(r) for r in results))

    def report(self, output_dict):
        """Print a table with the desired counts per assembly."""
        header = ('assembly', 'RS Created', 'RS Remapped', 'RS Deprecated', 'SS Clustered')
        output_table = [
            [asm, counts['clustered_variants_created'], counts.get('clustered_variants_remapped', 0),
             counts['clustered_variants_deprecated'], counts['submitted_variants_clustered']]
            for asm, counts in output_dict.items()
        ]
        # Remove all-zero rows from output
        output_table = [row for row in output_table if any(x > 0 for x in row[1:])]
        pretty_print(header, output_table)


if __name__ == '__main__':
    parser = ArgumentParser(description='')
    parser.add_argument('--settings-xml-file', required=True)
    parser.add_argument('--profile', required=True)
    parser.add_argument('--start_date', help='Start date in form YYYY-MM-DD, inclusive', required=True, type=str)
    parser.add_argument('--end_date', help='End date in form YYYY-MM-DD, exclusive', required=True, type=str)
    parser.add_argument('--job_tracker', help='Whether to use job tracker when gathering stats', action='store_true')

    args = parser.parse_args()
    counts = CountStats(args.profile, args.settings_xml_file)
    if args.job_tracker:
        counts.get_clustering_counts_between_dates_with_job_tracker(args.start_date, args.end_date)
    else:
        counts.get_clustering_counts_between_dates(args.start_date, args.end_date)
