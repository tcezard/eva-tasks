# Copyright 2022 EMBL - European Bioinformatics Institute
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

import argparse
import os
import re

import requests
from ebi_eva_common_pyutils.logger import logging_config

logger = logging_config.get_logger(__name__)


def prepare_processed_analyses_file(project, batch_size, processed_file_directory, target_file):
    total_analyses = total_analyses_in_project(project)
    logger.info(f"total analyses in project {project}: {total_analyses}")

    processed_analyses = get_processed_analyses_files(processed_file_directory)
    logger.info(f"no of processed analyses : {len(processed_analyses)}")
    processed_files_with_no_analyses = processed_analyses.copy()

    if not os.path.exists(target_file):
        logger.info(f"{target_file} does not exist and will be created")

    with open(target_file, 'a') as f:
        offset = 0
        limit = batch_size
        while offset < total_analyses:
            logger.info(f"Fetching ENA analyses from {offset} to  {offset + limit} (offset={offset}, limit={limit})")
            analyses_from_ena = get_analyses_from_ena(project, offset, limit)
            for analysis in analyses_from_ena:
                if analysis['run_ref'] in processed_analyses or analysis['analysis_accession'] in processed_analyses:
                    f.write(f"{analysis['analysis_accession']},{analysis['submitted_ftp']}\n")
                    processed_files_with_no_analyses.discard(analysis['run_ref'])

            offset = offset + limit

    print(f"Processed files for which no analyses were found in ENA : "
          f"total count = {len(processed_files_with_no_analyses)}, files = {processed_files_with_no_analyses}")


def total_analyses_in_project(project):
    count_url = f"https://www.ebi.ac.uk/ena/portal/api/filereportcount?accession={project}&result=analysis"
    response = requests.get(count_url)
    if response.status_code != 200:
        logger.error(f"Error fetching total analyses count for project {project}")
        response.raise_for_status()
    return response.json()


def get_processed_analyses_files(processed_file_directory):
    processed_analyses = set()
    for file in os.listdir(processed_file_directory):
        file_path = os.path.join(processed_file_directory, file)
        if os.path.islink(file_path):
            orig_path = os.readlink(file_path)
            logger.info(f"{file_path} is a symlink, reading files from original path {orig_path}")
            processed_analyses.update(get_processed_analyses_files(orig_path))
        elif os.path.isdir(file_path):
            logger.info(f"{file_path} is a directory, reading files from sub-directories")
            processed_analyses.update(get_processed_analyses_files(file_path))
        else:
            if (file.endswith("_filtered_vcf.gz") or file.endswith(".vcf") or file.endswith(".vcf.gz")) \
                    and not file.startswith("concat"):
                file_name = re.split('_|\.', file)[0]
                processed_analyses.add(file_name)

    return processed_analyses


def get_analyses_from_ena(project, offset, limit):
    analyses_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?result=analysis&accession={project}&offset={offset}" \
                   f"&limit={limit}&format=json&fields=run_ref,analysis_accession,submitted_ftp"
    response = requests.get(analyses_url)
    if response.status_code != 200:
        logger.error(f"Error fetching analyses info from ENA for {project}")
        response.raise_for_status()
    return response.json()


def main():
    parser = argparse.ArgumentParser(description='Download analyses for processing from the Covid-19 DP project',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--project", default='PRJEB45554', required=False,
                        help="project from which analyses needs to be downloaded")
    parser.add_argument("--batch-size", default=100000, required=False, help="batch size of ENA analyses download")
    parser.add_argument("--processed-file-directory", required=True,
                        help="full path to the directory where all the processed files are present")
    parser.add_argument("--target-file", required=True, help="full path to the target file that will be created")

    args = parser.parse_args()
    logging_config.add_stdout_handler()

    prepare_processed_analyses_file(args.project, args.batch_size, args.processed_file_directory, args.target_file)


if __name__ == "__main__":
    main()
