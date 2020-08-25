#!/usr/bin/env python
import csv
import glob
import logging
import os
import subprocess
import sys
from argparse import ArgumentParser
from collections import defaultdict

assigned_number_variant = 0
eva_accession_path = ''

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)-15s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def run_command_with_output(command_description, command, return_process_output=False):
    process_output = ""

    logger.info("Starting process: " + command_description)
    logger.info("Running command: " + command)

    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, universal_newlines=True,
                          shell=True) as process:
        for line in iter(process.stdout.readline, ''):
            line = str(line).rstrip()
            logger.info(line)
            if return_process_output:
                process_output += line + "\n"
        for line in iter(process.stderr.readline, ''):
            line = str(line).rstrip()
            logger.error(line)
    if process.returncode != 0:
        logger.error(command_description + " failed! Refer to the error messages for details.")
        raise subprocess.CalledProcessError(process.returncode, process.args)
    else:
        logger.info(command_description + " - completed successfully")
    if return_process_output:
        return process_output


def parse_input(input_file):
    """
    Parse the input with a dict reader and returns tthe data aggregated per taxid
    """
    data_per_taxid = defaultdict(list)
    data_per_taxid_and_assembly = defaultdict(list)

    with open(input_file) as open_file:
        reader = csv.DictReader(open_file, delimiter='\t')
        for row in reader:
            if row['Decision to include'] == 'yes':
                data_per_taxid_and_assembly[(row['Taxid'], row['Assembly'])].append(row)
                data_per_taxid[row['Taxid']].append(row)
    return data_per_taxid, data_per_taxid_and_assembly


def assign_tempmongo_host_round_robin(number_variants, total_number_of_variant):
    global assigned_number_variant
    assigned_number_variant += number_variants
    instance_number = int((assigned_number_variant // ((total_number_of_variant + 1) / 10)) + 1)
    return 'tempmongo-' + str(instance_number)


def retrieve_assembly(scientific_name, assembly_accession, assembly_dir, download_dir):
    fasta, report = None, None
    accession_dir = os.path.join(
        assembly_dir,
        scientific_name.lower().replace(' ', '_'),
        assembly_accession
    )

    fastas = glob.glob(os.path.join(accession_dir, '*.fa'))
    nucleotide_fastas = glob.glob(os.path.join(accession_dir, '*.fna'))
    reports = glob.glob(os.path.join(accession_dir, '*_assembly_report.txt'))
    if fastas and len(fastas) == 1:
        fasta = fastas[0]
    elif nucleotide_fastas and len(nucleotide_fastas) == 1:
        fasta = nucleotide_fastas[0]
    if reports and len(reports) == 1:
        report = reports[0]

    if fasta is None or report is None:
        report, fasta = download_assembly(scientific_name, assembly_accession, download_dir)
    return fasta, report


def download_assembly(scientific_name, assembly_accession, download_dir):
    python = os.path.join(eva_accession_path, "python-import-automation/bin/python")
    genome_downloader = os.path.join(eva_accession_path, "eva-accession/eva-accession-import-automation/genome_downloader.py")
    private_json = os.path.join(eva_accession_path, "private-config.json")

    output_dir = os.path.join(download_dir, scientific_name.lower().replace(' ', '_'))
    os.makedirs(output_dir, exist_ok=True)
    command = "{} {} -p {} -a {} -o {}""".format(python, genome_downloader, private_json, assembly_accession, output_dir)
    run_command_with_output('Retrieve assembly fasta and assembly report', command)

    assembly_report = glob.glob(os.path.join(output_dir, assembly_accession, "*_assembly_report.txt"))[0]
    assembly_fasta = os.path.join(output_dir, assembly_accession, assembly_accession + ".fa")
    return assembly_fasta, assembly_report


def count_variants(rows, only_variant_to_process=True):
    """Sum variants in the rows provided."""
    stay_unchanged = {row['Will stay unchanged'] for row in rows}
    return sum([
        int(row['Number Of Variants (submitted variants)'].replace(',', ''))
        for row in rows
        if not only_variant_to_process or ('no' in stay_unchanged and row['Source'] != 'DBSNP - filesystem')
    ])


def aggregate_list_of_species(input_file, assembly_dir, download_dir, output_assemblies_tsv, output_taxonmomy_tsv):
    """
    """
    data_per_taxid, data_per_taxid_and_assembly = parse_input(input_file)
    only_unmapped_tax_id = []
    unchanged_taxid = []
    dbsnps_only = []
    eva_only = []
    mixture = []
    taxid_to_tempmongo = {}
    total_number_variant = 0
    for taxid in data_per_taxid:
        rows = data_per_taxid[taxid]
        sources = {row['Source'] for row in rows}
        stay_unchanged = {row['Will stay unchanged'] for row in rows}

        number_variants = count_variants(rows)
        total_number_variant += number_variants

        source = None
        if len(sources) == 1:
            source = sources.pop()
        if source == 'DBSNP - filesystem':
            only_unmapped_tax_id.append(taxid)
        elif 'no' not in stay_unchanged:
            unchanged_taxid.append(taxid)
        elif source == 'EVA':
            eva_only.append(taxid)
        elif source == 'DBSNP':
            dbsnps_only.append(taxid)
        else:
            mixture.append(taxid)

    # Iterate a second time to assign a temp mongodb per taxonomy
    with open(output_taxonmomy_tsv, 'w') as open_output:
        headers = ['taxonomy_id', 'scientific_name', 'tempmongo_instance', 'should_be_copied',
                             'number_variants_to_process', 'total_num_variants']
        print('\t'.join([str(o) for o in headers]), file=open_output)

        for taxid in data_per_taxid:
            rows = data_per_taxid[taxid]
            scientific_name = {row['Scientific Name From Taxid'] for row in rows}.pop()
            number_variants_to_process = count_variants(rows)
            taxid_to_tempmongo[taxid] = assign_tempmongo_host_round_robin(number_variants_to_process,
                                                                          total_number_variant)

            if taxid in unchanged_taxid + only_unmapped_tax_id:
                to_copy = 'no'
                temp_mongo = ''
            else:
                to_copy = 'yes'
                temp_mongo = taxid_to_tempmongo[taxid]

            print('\t'.join([
                str(taxid),
                scientific_name,
                temp_mongo,
                to_copy,
                str(number_variants_to_process),
                str(count_variants(rows, False))
            ]), file=open_output)

    with open(output_assemblies_tsv, 'w') as open_output:
        headers = ['taxonomy_id', 'scientific_name', 'assembly', 'sources', 'fasta_path', 'report_path',
                   'tempmongo_instance', 'should_be_process', 'number_variants_to_process', 'total_num_variants']
        print('\t'.join([str(o) for o in headers]), file=open_output)
        for taxid, assembly in data_per_taxid_and_assembly:
            rows = data_per_taxid_and_assembly[(taxid, assembly)]
            scientific_name = {row['Scientific Name From Taxid'] for row in rows}.pop()
            number_variants = sum([int(row['Number Of Variants (submitted variants)'].replace(',', '')) for row in rows])

            if taxid in unchanged_taxid + only_unmapped_tax_id or assembly == 'Unmapped':
                to_process = 'no'
                fasta, report, temp_mongo = '', '', ''
            else:
                to_process = 'yes'
                fasta, report = retrieve_assembly(scientific_name, assembly, assembly_dir, download_dir)
                temp_mongo = taxid_to_tempmongo[taxid]

            out = [
                taxid,
                scientific_name,
                assembly,
                ','.join({row['Source'] for row in rows}),
                fasta,
                report,
                temp_mongo,
                to_process,
                number_variants,
                str(count_variants(rows, False))
            ]
            print('\t'.join([str(o) for o in out]), file=open_output)


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--input', help='Path to the file containing the taxonomies and assemblies', required=True)
    argparse.add_argument('--assembly_dir', help='Path to the directory containing pre-downloaded species assemblies', required=True)
    argparse.add_argument('--download_dir', help='Path to the directory where additional containing species assemblies will be downloaded', required=True)
    argparse.add_argument('--output_assemblies_tsv', help='Path to the tsv file that will contain the list of assemblies to process', required=True)
    argparse.add_argument('--output_taxonomy_tsv', help='Path to the tsv file that will contain the list of species to process', required=True)
    argparse.add_argument('--eva_accession_path', help='path to the directory that contain eva-accession code and private json file.')

    args = argparse.parse_args()

    global eva_accession_path
    if args.eva_accession_path:
        eva_accession_path = args.eva_accession_path
    
    aggregate_list_of_species(args.input, args.assembly_dir, args.download_dir, args.output_assemblies_tsv, args.output_taxonomy_tsv)


if __name__ == "__main__":
    main()
