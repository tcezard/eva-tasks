from argparse import ArgumentParser
from collections import defaultdict, Counter

from ebi_eva_common_pyutils.logger import AppLogger

from ebi_eva_common_pyutils.common_utils import pretty_print
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query


class CountStats(AppLogger):

    def __init__(self, profile, settings_file):
        self.profile = profile
        self.settings_file = settings_file

    def get_clustering_counts_between_dates(self, start_date, end_date):
        # Dates are strings of the form YYYY-MM-DD, to make my life easier
        # Inclusive of start, exclusive of end
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
                # TODO what about splits? creates a new RS but does not cluster any new SS
                if metric != 'clustered_variants_created':
                    results_by_assembly[assembly][metric] += count
                    continue
                if metric == 'clustered_variants_created' and count > 0:
                    # If submitted variants are newly clustered in this step, probably a true creation of RS
                    if v['submitted_variants_clustered'] > 0:
                        results_by_assembly[assembly][metric] += count
                    # If clustered variants are created and deprecated in a single step, probably a remediation
                    # These can be ignored  --> TODO confirm this
                    elif v['clustered_variants_deprecated'] > 0:
                        continue
                    # Otherwise this is probably creation of clustered variants due to remapping
                    else:
                        # note for remapped this assembly is the target
                        results_by_assembly[assembly]['clustered_variants_remapped'] += count

        self.report(results_by_assembly)

    def report(self, output_dict):
        """Print a table with the desired counts per assembly."""
        header = ('assembly', 'RS Created', 'RS Remapped', 'RS Deprecated', 'SS Clustered')
        output_table = [
            [asm, counts['clustered_variants_created'], counts['clustered_variants_remapped'],
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

    args = parser.parse_args()
    counts = CountStats(args.profile, args.settings_xml_file)
    counts.get_clustering_counts_between_dates(args.start_date, args.end_date)
