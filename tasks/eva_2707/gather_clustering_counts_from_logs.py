import argparse
import csv
import glob
import os
from collections import defaultdict
from csv import excel_tab
from datetime import datetime

from ebi_eva_common_pyutils.logger import logging_config

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def gather_count_from_logs(clustering_dir, output_file):
    pass
    # Assume the directory structure:
    # clustering_dir --> <scientific_name_taxonomy_id> --> <assembly_accession> --> cluster_<date>.log

    all_log_pattern = os.path.join(clustering_dir, '*', 'GCA*', 'cluster_*.log')
    all_log_files = glob.glob(all_log_pattern)
    metrics_per_species = defaultdict(dict)
    for log_file in all_log_files:
        logger.info('Parse log file: ' + log_file)
        scientific_name, taxid, assembly_accession, log_date = parse_log_file_path(log_file)
        result_dict = parse_one_log(log_file)
        truncated = detect_missing_data(result_dict)
        if truncated:
            metrics_per_species[taxid]['truncated'] = 'Yes'
        relevant_metrics = extract_relevant_metrics_from_one_log(result_dict)
        metrics_per_species[taxid]['taxid'] = taxid
        metrics_per_species[taxid]['scientific_name'] = scientific_name
        metrics_per_species[taxid]['assembly_accession'] = assembly_accession
        if 'start_date' not in metrics_per_species[taxid] or \
                log_date < metrics_per_species[taxid]['start_date']:
            metrics_per_species[taxid]['start_date'] = log_date
        for metric in relevant_metrics:
            if metric in metrics_per_species[taxid]:
                metrics_per_species[taxid][metric] += int(relevant_metrics[metric])
            else:
                metrics_per_species[taxid][metric] = int(relevant_metrics[metric])

    with open(output_file, 'w') as open_file:
        fieldnames = ['taxid', 'scientific_name', 'assembly_accession', 'start_date', 'clustered_variant_remapped',
                      'clustered_variants_created', 'submitted_variants_clustered', 'merge_operations',
                      'split_operations', 'truncated']
        writer = csv.DictWriter(open_file, fieldnames=fieldnames, dialect=excel_tab)
        writer.writeheader()
        for taxid in metrics_per_species:
            writer.writerow(metrics_per_species[taxid])


def detect_missing_data(result_dict):
    truncated = False
    for step in steps:
        if step not in result_dict:
            logger.warning('Missing ' + step + ' data ')
            result_dict[step] = {}
        else:
            if list(result_dict[step].keys()) == ['item_written']:
                # Missing all the count but the step has started so we could have lost some information.
                logger.error('Truncated data in step ' + step)
                truncated = True
    return truncated


def parse_log_file_path(log_file_path):
    scientific_name_taxonomy_id, assembly_accession, file_name = log_file_path.split('/')[-3:]
    scientific_name = '_'.join(scientific_name_taxonomy_id.split('_')[:-1])
    taxid = scientific_name_taxonomy_id.split('_')[-1]
    date = datetime.strptime(file_name.split('.')[0].split('_')[-1], '%Y%m%d%H%M%S')  # 20220112170519
    return scientific_name, taxid, assembly_accession, date


def parse_one_log(log_file):
    # identify the clustering step
    # identify the start end of run
    # find the count lines and extract metrics
    results = {}
    with open(log_file) as open_file:
        current_step = None
        for line in open_file:
            sp_line = line.strip().split()
            if len(sp_line) < 8:
                continue
            if sp_line[7] == 'u.a.e.e.a.c.b.l.GenericProgressListener':
                current_step = sp_line[9].rstrip(':')
                if current_step not in results:
                    results[current_step] = {}
                if len(sp_line) > 17:
                    results[current_step]['item_written'] = sp_line[17]
            elif sp_line[7] == 'u.a.e.eva.metrics.metric.MetricCompute' and sp_line[9] == 'Count{id=null,':
                metric = sp_line[12].split("'")[1]
                count = sp_line[13].split('=')[1].rstrip('}')
                results[current_step][metric] = count
    return results

steps = [
        'CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP', 'PROCESS_RS_MERGE_CANDIDATES_STEP',
        'PROCESS_RS_SPLIT_CANDIDATES_STEP', 'CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP'
    ]

metrics = [
    'clustered_variants_created',
    'clustered_variants_updated',
    'clustered_variants_rs_split',
    'clustered_variants_merged',
    'clustered_variants_merge_operations',
    'submitted_variants_clustered',
    'submitted_variants_kept_unclustered',
    'submitted_variants_updated_rs',
    'submitted_variants_update_operations',
]


def extract_relevant_metrics_from_one_log(parsed_log_dict):
    """Select only the metrics that will matter and rename them"""
    # Steps CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP and PROCESS_RS_SPLIT_CANDIDATES_STEP are run after
    # PROCESS_RS_SPLIT_CANDIDATES_STEP and the metrics are not reset between the steps
    for metric in metrics:

        parsed_log_dict['CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP'][metric] = max(
            int(parsed_log_dict.get('CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP', {}).get(metric, 0)) -
            int(parsed_log_dict.get('PROCESS_RS_SPLIT_CANDIDATES_STEP', {}).get(metric, 0)),
            0
        )

        parsed_log_dict['PROCESS_RS_MERGE_CANDIDATES_STEP'][metric] = max(
            int(parsed_log_dict.get('PROCESS_RS_SPLIT_CANDIDATES_STEP', {}).get(metric, 0)) -
            int(parsed_log_dict.get('PROCESS_RS_MERGE_CANDIDATES_STEP', {}).get(metric, 0)),
            0
        )


    # Only metrics that matters for the first step
    clustered_variant_remapped = parsed_log_dict.get('CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP', {}).get('clustered_variants_created', 0)
    merge_operations = parsed_log_dict.get('PROCESS_RS_MERGE_CANDIDATES_STEP', {}).get('clustered_variants_merge_operations', 0)
    split_operations = parsed_log_dict.get('PROCESS_RS_SPLIT_CANDIDATES_STEP', {}).get('clustered_variants_merge_operations', 0)
    clustered_variants_created = parsed_log_dict.get('CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP', {}).get('clustered_variants_created', 0)
    submitted_variants_clustered = parsed_log_dict.get('CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP', {}).get('submitted_variants_clustered', 0)
    return {
        'clustered_variant_remapped': clustered_variant_remapped,
        'merge_operations': merge_operations,
        'split_operations': split_operations,
        'clustered_variants_created': clustered_variants_created,
        'submitted_variants_clustered': submitted_variants_clustered
    }


def main():
    parser = argparse.ArgumentParser(
        description='Parse all the clustering logs and report number of event created in each one')
    parser.add_argument("--clustering_root_path", type=str,
                        help="base directory where all the clustering was run.", required=True)
    parser.add_argument("--output_csv", type=str,
                        help="path to the output .", required=True)

    args = parser.parse_args()
    gather_count_from_logs(args.clustering_root_path, args.output_csv)


if __name__ == '__main__':
    main()
