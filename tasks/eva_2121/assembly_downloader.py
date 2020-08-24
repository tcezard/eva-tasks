# Copyright 2019 EMBL - European Bioinformatics Institute
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
import logging
import os
from ftplib import FTP
import re
from urllib.parse import urlparse

import wget
import pandas
import sys
import urllib.request
import re
import json
from retry import retry


logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)-15s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def get_assembly_report_url(assembly_accession):
    if re.match(r"^GC[F|A]_\d+\.\d+$", assembly_accession) is None:
        raise Exception('Invalid assembly accession: it has to be in the form of '
                        'GCF_XXXXXXXXX.X or GCA_XXXXXXXXX.X where X is a number')

    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=600)
    ftp.login()

    genome_folder = 'genomes/all/' + '/'.join([assembly_accession[0:3], assembly_accession[4:7],
                                               assembly_accession[7:10], assembly_accession[10:13]]) + '/'
    ftp.cwd(genome_folder)

    all_genome_subfolders = []
    ftp.retrlines('NLST', lambda line: all_genome_subfolders.append(line))

    genome_subfolders = [folder for folder in all_genome_subfolders if assembly_accession in folder]

    if len(genome_subfolders) != 1:
        raise Exception('more than one folder matches the assembly accession: ' + str(genome_subfolders))

    ftp.cwd(genome_subfolders[0])
    genome_files = []
    ftp.retrlines('NLST', lambda line: genome_files.append(line))
    ftp.quit()

    assembly_reports = [genome_file for genome_file in genome_files if 'assembly_report' in genome_file]
    if len(assembly_reports) != 1:
        raise Exception('more than one file has "assembly_report" in its name: ' + str(assembly_reports))

    return 'ftp://' + 'ftp.ncbi.nlm.nih.gov' + '/' + genome_folder + genome_subfolders[0] + '/' + assembly_reports[0]


def download_assembly_report(detination_dir, assembly_accession):
    print('Searching assembly report for {} ...'.format(assembly_accession))
    url = get_assembly_report_url(assembly_accession)
    asm_report_file = os.path.join(detination_dir, os.path.basename(urlparse(url).path))
    downloaded_already = os.path.isfile(asm_report_file)
    if downloaded_already:
        print('Assembly report was already downloaded, skipping download ...')
    else:
        print('Downloading assembly report ...')
        download_file_from_ftp(url, detination_dir)
    return asm_report_file


def download_file_from_ftp(url, destination_dir):
    return wget.download(url, out=destination_dir)


def build_fasta_from_assembly_report(assembly_report_path, assembly_accession, eutils_api_key, directory_path, clear):
    fasta_path = directory_path + assembly_accession + '.fa'
    if clear:
        os.remove(fasta_path)
    written_contigs = get_written_contigs(fasta_path)
    for chunk in pandas.read_csv(assembly_report_path, skiprows=get_header_line_index(assembly_report_path),
                                 dtype=str, sep='\t', chunksize=100):
        new_contigs_in_chunk = process_chunk(chunk, written_contigs, eutils_api_key, fasta_path)
        written_contigs.extend(new_contigs_in_chunk)


def get_written_contigs(fasta_path):
    try:
        written_contigs = []
        match = re.compile(r'>(.*?)\s')
        with open(fasta_path, 'r') as file:
            for line in file:
                written_contigs.extend(match.findall(line))
        return written_contigs
    except FileNotFoundError:
        logger.info('FASTA file does not exists, starting from scratch')
        return []


def get_header_line_index(assembly_report_path):
    try:
        with open(assembly_report_path) as assembly_report_file_handle:
            return next(line_index for line_index, line in enumerate(assembly_report_file_handle)
                        if line.lower().startswith("# sequence-name") and "sequence-role" in line.lower())
    except StopIteration:
        raise Exception("Could not locate header row in the assembly report!")


def process_chunk(assembly_report_dataframe, written_contigs, eutils_api_key, fasta_path):
    new_contigs = []
    for index, row in assembly_report_dataframe.iterrows():
        genbank_accession = row['GenBank-Accn']
        refseq_accession = row['RefSeq-Accn']
        relationship = row['Relationship']
        accession = genbank_accession
        if relationship != '=' and genbank_accession == 'na':
            accession = refseq_accession
        if written_contigs is not None and (genbank_accession in written_contigs or refseq_accession in written_contigs):
            logger.info('Accession ' + accession + ' already in the FASTA file, don\'t need to be downloaded')
            continue
        new_contig = get_sequence_from_ncbi(accession, fasta_path, eutils_api_key)
        if new_contig is not None:
            new_contigs.append(new_contig)
    return new_contigs


@retry(tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def get_sequence_from_ncbi(accession, fasta_path, eutils_api_key):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + accession + \
           '&rettype=fasta&retmode=text&api_key=' + eutils_api_key + '&tool=eva&email=eva-dev@ebi.ac.uk'
    logger.info('Downloading ' + accession)
    sequence_tmp_path = os.path.dirname(fasta_path) + '/' + accession
    urllib.request.urlretrieve(url, sequence_tmp_path)
    with open(sequence_tmp_path) as sequence:
        first_line = sequence.readline()
        second_line = sequence.readline()
        if not second_line:
            logger.info('FASTA sequence not available for ' + accession)
            os.remove(sequence_tmp_path)
        else:
            contatenate_sequence_to_fasta(fasta_path, sequence_tmp_path)
            logger.info(accession + " downloaded and added to FASTA sequence")
            os.remove(sequence_tmp_path)
            return accession


def contatenate_sequence_to_fasta(fasta_path, sequence_path):
    with open(fasta_path, 'a+') as fasta:
        with open(sequence_path) as sequence:
            for line in sequence:
                # Check that the line is not empty
                if line.strip():
                    fasta.write(line)


@retry(tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def download_assembly_report(assembly_accession, directory_path):
    assembly_report_url = get_assembly_report_url(assembly_accession)
    assembly_report_path = directory_path + os.path.basename(assembly_report_url)
    os.makedirs(os.path.dirname(assembly_report_path), exist_ok=True)
    urllib.request.urlretrieve(assembly_report_url, assembly_report_path)
    urllib.request.urlcleanup()
    return assembly_report_path


def build_output_directory_path(assembly_accession, private_config_args, output_directory):
    if output_directory is not None:
        directory_path = output_directory
    else:
        eva_root_dir = private_config_args['eva_root_dir']
        directory_path = eva_root_dir
    output_directory += os.path.sep + assembly_accession + os.path.sep
    logger.info('Files will be downloaded in ' + directory_path)
    return directory_path


def download_assembly(assembly_accession, species_name, output_directory, eutils_api_key, clear):
    directory_path = os.path.join(output_directory, species_name, assembly_accession)
    assembly_report_path = download_assembly_report(assembly_accession, directory_path)
    build_fasta_from_assembly_report(assembly_report_path, assembly_accession, eutils_api_key, directory_path, clear)

