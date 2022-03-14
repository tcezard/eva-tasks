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

import requests
from ebi_eva_common_pyutils.logger import logging_config

logger = logging_config.get_logger(__name__)


def prepare_processed_analysis_file(project, batch_size, processed_file_directory, target_file):
    total_analysis = total_analysis_in_project(project)
    logger.info(f"total analysis in project {project}: {total_analysis}")

    processed_analysis = get_processed_analysis_files(processed_file_directory)

    if not os.path.exists(target_file):
        logger.info(f"{target_file} does not exist and will be created")

    with open(target_file, 'a') as f:
        offset = 0
        limit = batch_size
        while offset < total_analysis:
            analysis_from_ena = get_analysis_from_ena(project, offset, limit)
            for analysis in analysis_from_ena:
                if analysis['run_ref'] in processed_analysis:
                    logger.info(f"Analysis file {analysis['run_ref']} found in already processed list")
                    f.write(f"{analysis['analysis_accession']},{analysis['submitted_ftp']}\n")
            offset = offset + limit


def total_analysis_in_project(project):
    count_url = f"https://www.ebi.ac.uk/ena/portal/api/filereportcount?accession={project}&result=analysis"
    response = requests.get(count_url)
    if response.status_code != 200:
        logger.error(f"Error fetching total analysis count for project {project}")
        response.raise_for_status()
    return response.json()


def get_processed_analysis_files(processed_file_directory):
    processed_analysis = set()
    for file in os.listdir(processed_file_directory):
        file_path = os.path.join(processed_file_directory, file)
        if os.path.islink(file_path):
            orig_path = os.readlink(file_path)
            logger.info(f"{file_path} is a symlink, reading files from original path {orig_path}")
            processed_analysis.update(get_processed_analysis_files(orig_path))
        elif os.path.isdir(file_path):
            logger.info(f"{file_path} is a directory, reading files from sub-directories")
            processed_analysis.update(get_processed_analysis_files(file_path))
        else:
            if "_filtered_vcf.gz" in file:
                file_name = file.split("_")[0]
                processed_analysis.add(file_name)

    return processed_analysis


def get_analysis_from_ena(project, offset, limit):
    analysis_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?result=analysis&accession={project}&offset={offset}" \
                   f"&limit={limit}&format=json&fields=run_ref,analysis_accession,submitted_ftp"
    response = requests.get(analysis_url)
    if response.status_code != 200:
        logger.error(f"Error fetching analysis info from ENA for {project}")
        response.raise_for_status()
    return response.json()


def main():
    parser = argparse.ArgumentParser(description='Download analysis for processing from the Covid-19 DP project',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--project", default='PRJEB45554', required=False,
                        help="project from which analysis needs to be downloaded")
    parser.add_argument("--batch-size", default=10000, required=False, help="batch size of ENA analysis download")
    parser.add_argument("--processed-file-directory", required=True,
                        help="full path to the directory where all the processed files are present")
    parser.add_argument("--target-file", required=True, help="full path to the target file that will be created")

    args = parser.parse_args()
    logging_config.add_stdout_handler()

    prepare_processed_analysis_file(args.project, args.batch_size, args.processed_file_directory, args.target_file)


if __name__ == "__main__":
    main()