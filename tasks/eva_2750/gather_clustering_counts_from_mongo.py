import argparse
import csv
import glob
import os
from collections import defaultdict
from csv import excel_tab
from datetime import datetime

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase


logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def gather_count_from_mongo(clustering_dir, output_file, mongo_source):
    # Assume the directory structure:
    # clustering_dir --> <scientific_name_taxonomy_id> --> <assembly_accession> --> cluster_<date>.log_dict

    all_log_pattern = os.path.join(clustering_dir, '*', 'GCA*', 'cluster_*.log')
    all_log_files = glob.glob(all_log_pattern)
    ranges_per_assembly = defaultdict(dict)
    for log_file in all_log_files:
        logger.info('Parse log_dict file: ' + log_file)
        scientific_name, taxid, assembly_accession, log_date = parse_log_file_path(log_file)
        log_metric_date_range = parse_one_log(log_file)

        if assembly_accession not in ranges_per_assembly:
            ranges_per_assembly[assembly_accession] = defaultdict(dict)

        # new_remapped_current_rs
        if 'CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP' in log_metric_date_range:
            ranges_per_assembly[assembly_accession]['new_remapped_current_rs'][log_file] = {
                'from': log_metric_date_range['CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP'],
                'to': log_metric_date_range['CLEAR_RS_MERGE_AND_SPLIT_CANDIDATES_STEP']
                if "CLEAR_RS_MERGE_AND_SPLIT_CANDIDATES_STEP" in log_metric_date_range
                else log_metric_date_range["last_timestamp"]
            }
        # new_current_rs
        if 'CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP' in log_metric_date_range:
            ranges_per_assembly[assembly_accession]['new_current_rs'][log_file] = {
                'from': log_metric_date_range['CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP'],
                'to': log_metric_date_range['CLUSTER_UNCLUSTERED_VARIANTS_JOB']["completed"]
                if "completed" in log_metric_date_range['CLUSTER_UNCLUSTERED_VARIANTS_JOB']
                else log_metric_date_range["last_timestamp"]
            }
        # merged_rs
        if 'PROCESS_RS_MERGE_CANDIDATES_STEP' in log_metric_date_range:
            ranges_per_assembly[assembly_accession]['merged_rs'][log_file] = {
                'from': log_metric_date_range['PROCESS_RS_MERGE_CANDIDATES_STEP'],
                'to': log_metric_date_range['PROCESS_RS_SPLIT_CANDIDATES_STEP']
                if "PROCESS_RS_SPLIT_CANDIDATES_STEP" in log_metric_date_range
                else log_metric_date_range["last_timestamp"]
            }
        # split_rs
        if 'PROCESS_RS_SPLIT_CANDIDATES_STEP' in log_metric_date_range:
            ranges_per_assembly[assembly_accession]['split_rs'][log_file] = {
                'from': log_metric_date_range['PROCESS_RS_SPLIT_CANDIDATES_STEP'],
                'to': log_metric_date_range['CLEAR_RS_MERGE_AND_SPLIT_CANDIDATES_STEP']
                if "CLEAR_RS_MERGE_AND_SPLIT_CANDIDATES_STEP" in log_metric_date_range
                else log_metric_date_range["last_timestamp"]
            }
        # new_ss_clustered
        if 'CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP' in log_metric_date_range:
            ranges_per_assembly[assembly_accession]['new_ss_clustered'][log_file] = {
                'from': log_metric_date_range['CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP'],
                'to': log_metric_date_range['CLUSTER_UNCLUSTERED_VARIANTS_JOB']["completed"]
                if "completed" in log_metric_date_range['CLUSTER_UNCLUSTERED_VARIANTS_JOB']
                else log_metric_date_range["last_timestamp"]
            }

    metrics_per_assembly = defaultdict(dict)
    for asm, count_dict in ranges_per_assembly.items():
        new_remapped_current_rs, new_current_rs, merged_rs, split_rs, new_ss_clustered = 0, 0, 0, 0, 0
        for metric, log_dict in count_dict.items():
            expressions = []
            for log_name, query_range in log_dict.items():
                expressions.append({"createdDate": {"$gt": query_range["from"], "$lt": query_range["to"]}})

            date_range_filter = expressions
            if metric == 'new_remapped_current_rs':
                filter_criteria = {'asm': asm, '$or': date_range_filter}
                new_remapped_current_rs = query_mongo(mongo_source, filter_criteria, metric)
                logger.info(f'{metric} = {new_remapped_current_rs})')
            elif metric == 'new_current_rs':
                filter_criteria = {'asm': asm, '$or': date_range_filter}
                new_current_rs = query_mongo(mongo_source, filter_criteria, metric)
                logger.info(f'{metric} = {new_current_rs})')
            elif metric == 'merged_rs':
                filter_criteria = {'inactiveObjects.asm': asm, 'eventType': 'MERGED',
                                   '$or': date_range_filter}
                merged_rs = query_mongo(mongo_source, filter_criteria, metric)
                logger.info(f'{metric} = {merged_rs})')
            elif metric == 'split_rs':
                filter_criteria = {'inactiveObjects.asm': asm, 'eventType': 'RS_SPLIT',
                                   '$or': date_range_filter}
                split_rs = query_mongo(mongo_source, filter_criteria, metric)
                logger.info(f'{metric} = {split_rs})')
            elif metric == 'new_ss_clustered':
                filter_criteria = {'inactiveObjects.seq': asm, 'eventType': 'UPDATED',
                                   '$or': date_range_filter}
                new_ss_clustered = query_mongo(mongo_source, filter_criteria, metric)
                logger.info(f'{metric} = {new_ss_clustered})')

        metrics_per_assembly[asm]["assembly_accession"] = asm
        metrics_per_assembly[asm]["new_remapped_current_rs"] = new_remapped_current_rs
        metrics_per_assembly[asm]["new_current_rs"] = new_current_rs
        metrics_per_assembly[asm]["merged_rs"] = merged_rs
        metrics_per_assembly[asm]["split_rs"] = split_rs
        metrics_per_assembly[asm]["new_ss_clustered"] = new_ss_clustered

    with open(output_file, 'w') as open_file:
        fieldnames = ['assembly_accession', 'new_remapped_current_rs', 'new_current_rs', 'merged_rs', 'split_rs',
                      'new_ss_clustered']
        writer = csv.DictWriter(open_file, fieldnames=fieldnames, dialect=excel_tab)
        writer.writeheader()
        for asm in metrics_per_assembly:
            writer.writerow(metrics_per_assembly[asm])


def query_mongo(mongo_source, filter_criteria, metric):
    total_count = 0
    for collection_name in collections[metric]:
        logger.info(f'Querying mongo: db.{collection_name}.countDocuments({filter_criteria})')
        collection = mongo_source.mongo_handle[mongo_source.db_name][collection_name]
        count = collection.count_documents(filter_criteria)
        total_count += count
        logger.info(f'{count}')
    return total_count


collections = {
    "new_remapped_current_rs": [
        "clusteredVariantEntity",
        "dbsnpClusteredVariantEntity"
    ],
    "new_current_rs": [
        "clusteredVariantEntity"
    ],
    "merged_rs": [
        "clusteredVariantOperationEntity",
        "dbsnpClusteredVariantOperationEntity"
    ],
    "split_rs": [
        "clusteredVariantOperationEntity",
        "dbsnpClusteredVariantOperationEntity"
    ],
    "new_ss_clustered": [
        "submittedVariantOperationEntity",
        "dbsnpSubmittedVariantOperationEntity"
    ]
}


def parse_log_file_path(log_file_path):
    scientific_name_taxonomy_id, assembly_accession, file_name = log_file_path.split('/')[-3:]
    scientific_name = '_'.join(scientific_name_taxonomy_id.split('_')[:-1])
    taxid = scientific_name_taxonomy_id.split('_')[-1]
    date = datetime.strptime(file_name.split('.')[0].split('_')[-1], '%Y%m%d%H%M%S')  # 20220112170519
    return scientific_name, taxid, assembly_accession, date


def parse_one_log(log_file):
    # identify the clustering job/step
    # identify the start end of run
    # find the count lines and extract metrics
    results = {}
    with open(log_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            if len(sp_line) < 8:
                continue

            # get last timestamp
            try:
                timestamp = datetime.strptime(f"{sp_line[0]}T{sp_line[1]}Z", '%Y-%m-%dT%H:%M:%S.%fZ')
            except ValueError:
                pass

            # Jobs
            if sp_line[7] == 'o.s.b.c.l.support.SimpleJobLauncher':
                if sp_line[12] == "launched" or sp_line[12] == "completed":
                    current_job = sp_line[11].rstrip(']').lstrip('[name=')
                    job_status = sp_line[12]
                    if current_job not in results:
                        results[current_job] = {}
                    if job_status not in results[current_job]:
                        results[current_job][job_status] = {}
                    results[current_job][job_status] = timestamp
            # Steps
            if sp_line[7] == 'o.s.batch.core.job.SimpleStepHandler':
                current_step = sp_line[11].rstrip(']').lstrip('[')
                if current_step not in results:
                    results[current_step] = {}
                results[current_step] = timestamp

    results["last_timestamp"] = timestamp
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Parse all the clustering logs to get date ranges and query mongo to get metrics counts')
    parser.add_argument("--clustering_root_path", type=str,
                        help="base directory where all the clustering was run.", required=True)
    parser.add_argument("--output_csv", type=str,
                        help="path to the output .", required=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)

    args = parser.parse_args()
    mongo_source = MongoDatabase(uri=args.mongo_source_uri, secrets_file=args.mongo_source_secrets_file,
                                 db_name="eva_accession_sharded")
    gather_count_from_mongo(args.clustering_root_path, args.output_csv, mongo_source)


if __name__ == '__main__':
    main()
