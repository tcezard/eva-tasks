#!/usr/bin/env python
import csv
import glob
import json
import logging
import os
import shutil
from argparse import ArgumentParser
from collections import defaultdict
from operator import itemgetter

import psycopg2
from ebi_eva_common_pyutils.assembly import NCBIAssembly
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query

assigned_number_variant = 0
number_of_tempmongo_instances = 10
eva_accession_path = ''

logger = log_cfg.get_logger(__name__)


def parse_application_properties(app_prop):
    properties_dict = { }
    with open(app_prop) as open_file:
        for line in open_file:
            if not line.strip() or line.strip().startswith('#'):
                continue
            sp_line = line.strip().split('=')
            properties_dict[sp_line[0].strip()] = sp_line[1].strip()
    return properties_dict


def find_and_parse_properties(properties_path):
    property_files = glob.glob(os.path.join(properties_path, '*', 'release*', '*', '*application.properties'))
    property_files.extend(glob.glob(os.path.join(properties_path, '*', 'release*', '*application.properties')))
    all_properties_files_per_accession = defaultdict(list)
    for property_file in property_files:
        properties_dict = parse_application_properties(property_file)
        # if os.path.isfile(properties_dict['parameters.assemblyReportUrl']):
        all_properties_files_per_accession[properties_dict['parameters.assemblyAccession']].append({
            'content': properties_dict,
            'date_modified': os.path.getmtime(property_file),
            'file_path': property_file
        })

    # now keep the most recent one
    most_recent_per_accession = {}
    for accession in all_properties_files_per_accession:
        logger.debug('Out of the following files:\n' + '\n'.join([f['file_path'] for f in all_properties_files_per_accession[accession]]))
        most_recent_per_accession[accession] = sorted(all_properties_files_per_accession[accession], key=itemgetter('date_modified'))[-1]
        logger.debug('Most recent is :' + most_recent_per_accession[accession]['file_path'])

    return most_recent_per_accession


def resolve_fasta_and_report_path_from_release1(properties_path):
    property_per_accession = find_and_parse_properties(properties_path)
    current_dir = os.getcwd()
    paths_per_accession = {}
    for accession in property_per_accession:
        os.chdir(os.path.dirname(property_per_accession[accession]['file_path']))
        fasta = os.path.abspath(property_per_accession[accession]['content']['parameters.fasta'])
        report = property_per_accession[accession]['content']['parameters.assemblyReportUrl']
        if report.startswith('file:'):
            report = os.path.abspath(report[5:])
        else:
            report = None
        paths_per_accession[accession] = (fasta, report)
    os.chdir(current_dir)
    return paths_per_accession


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
    instance_number = int(
        (assigned_number_variant // ((total_number_of_variant + 1) / number_of_tempmongo_instances)) + 1
    )
    return 'tempmongo-' + str(instance_number)


def look_for_assembly_file(assembly_path, scientific_name, assembly_accession, pattern):
    accession_dir = os.path.join(
        assembly_path,
        scientific_name.lower().replace(' ', '_'),
        assembly_accession
    )
    file_paths = glob.glob(os.path.join(accession_dir, pattern))
    if len(file_paths) == 1:
        logger.debug('Found a single file in %s with pattern %s: %s', accession_dir, pattern, file_paths[0])
        return file_paths[0]
    else:
        logger.debug('Cannot find a single file in %s with pattern %s. Found %s', accession_dir, pattern, len(file_paths))


def search_for_assembly_in_dir(scientific_name, assembly_accession, assembly_search_paths):
    fasta, report = None, None
    for assembly_search_path in assembly_search_paths:

        if not fasta:
            fasta = look_for_assembly_file(assembly_search_path, scientific_name, assembly_accession, assembly_accession + '_custom.fa')
        if not fasta:
            fasta = look_for_assembly_file(assembly_search_path, scientific_name, assembly_accession, assembly_accession + '.fa')
        if not fasta:
            fasta = look_for_assembly_file(assembly_search_path, scientific_name, assembly_accession, '*.fa')
        if not fasta:
            fasta = look_for_assembly_file(assembly_search_path, scientific_name, assembly_accession, '*.fna')
        if not report:
            report = look_for_assembly_file(assembly_search_path, scientific_name, 'dbsnp_import', '*_custom_assembly_report.txt')
        if not report:
            report = look_for_assembly_file(assembly_search_path, scientific_name, 'dbsnp_import', '*_assembly_report_CUSTOM.txt')
        if not report:
            report = look_for_assembly_file(assembly_search_path, scientific_name, assembly_accession, '*_assembly_report.txt')

    if fasta and report:
        logger.info('Found Assembly fasta and report for %s, %s', scientific_name, assembly_accession)
        logger.info('fasta file: %s', fasta)
        logger.info('report file: %s', report)

    return fasta, report


def download_assembly(scientific_name, assembly_accession, download_dir, assembly_report=None):
    private_json = os.path.join(eva_accession_path, "private-config.json")
    with open(private_json) as private_config_file_handle:
        config = json.load(private_config_file_handle)
        eutils_api_key = config['eutils_api_key']

    assembly = NCBIAssembly(assembly_accession, scientific_name, download_dir, eutils_api_key)
    if assembly_report:
        shutil.copyfile(assembly_report, assembly.assembly_report_path)
    assembly.download_or_construct()
    return assembly.assembly_fasta_path, assembly.assembly_report_path


def count_variants(rows, only_variant_to_process=True):
    """Sum variants in the rows provided."""
    stay_unchanged = {row['Will stay unchanged'] for row in rows}
    return sum([
        int(row['Number Of Variants (submitted variants)'].replace(',', ''))
        for row in rows
        if not only_variant_to_process or ('no' in stay_unchanged and row['Source'] != 'DBSNP - filesystem')
    ])


def prepare_location_of_report_and_fasta(path_per_assemblies_from_release1, assembly, scientific_name, download_dir,
                                         assembly_dirs, release2_reference_folder):
    # First check previous release for a suitable assembly
    fasta, report = path_per_assemblies_from_release1.get(assembly, (None, None))

    # If not found then search for existing assembly in different directories
    if not report or not os.path.isfile(report):
        fasta, report = search_for_assembly_in_dir(scientific_name, assembly, assembly_dirs)
    else:
        logger.info('Found Assembly fasta and report for %s, %s In RELEASE1', scientific_name, assembly)
        logger.info('fasta file: %s', fasta)
        logger.info('report file: %s', report)

    # if the report is not there then download the fasta ans the report
    if not report or not os.path.isfile(report):
        fasta, report = download_assembly(scientific_name, assembly, download_dir)
    # If the report is there then download the fasta based on the report
    elif not fasta or not os.path.isfile(fasta):
        fasta, report = download_assembly(scientific_name, assembly, download_dir, report)

    assembly_folder = os.path.join(release2_reference_folder, scientific_name.lower().replace(' ', '_'), assembly)
    os.makedirs(assembly_folder, exist_ok=True)
    shutil.copyfile(fasta, os.path.join(assembly_folder, assembly + '.fa'))
    shutil.copyfile(report, os.path.join(assembly_folder, assembly + '_assembly_report.txt'))
    return os.path.join(assembly_folder, assembly + '.fa'), os.path.join(assembly_folder, assembly + '_assembly_report.txt')


def get_dbsnp_database_name(pg_conn, scientific_name):
    query = "select database_name from dbsnp_ensembl_species.import_progress where scientific_name='{}';"
    database_names = list(set([database_name for database_name, in get_all_results_for_query(pg_conn, query.format(scientific_name))]))
    assert len(database_names) < 2, "Species {} has more than one database associated: {}".format(scientific_name, database_names)
    if database_names:
        return database_names[0]


def aggregate_list_of_species(input_file, properties_path, assembly_dirs, download_dir, release2_reference_folder,
                              output_assemblies_tsv, output_taxonmomy_tsv, private_config_xml_file):
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

    path_per_assemblies_from_release1 = resolve_fasta_and_report_path_from_release1(properties_path)

    pg_conn = psycopg2.connect(
        get_pg_metadata_uri_for_eva_profile("development", private_config_xml_file), user="evadev"
    )

    with open(output_assemblies_tsv, 'w') as open_output:
        headers = ['taxonomy_id', 'scientific_name', 'assembly', 'sources', 'fasta_path', 'report_path',
                   'tempmongo_instance', 'should_be_process', 'number_variants_to_process', 'total_num_variants',
                   'dbsnp_database_name', 'release_version']
        print('\t'.join([str(o) for o in headers]), file=open_output)
        for taxid, assembly in data_per_taxid_and_assembly:
            rows = data_per_taxid_and_assembly[(taxid, assembly)]
            scientific_name = {row['Scientific Name From Taxid'] for row in rows}.pop()
            dbsnps_database_name = get_dbsnp_database_name(pg_conn, scientific_name) or ''

            if taxid in unchanged_taxid + only_unmapped_tax_id or assembly == 'Unmapped':
                to_process = 'no'
                fasta, report, temp_mongo = '', '', ''
            else:
                to_process = 'yes'
                fasta, report = prepare_location_of_report_and_fasta(
                    path_per_assemblies_from_release1, assembly, scientific_name,  download_dir,
                    assembly_dirs, release2_reference_folder
                )
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
                str(count_variants(rows, True)),
                str(count_variants(rows, False)),
                '2',
                dbsnps_database_name
            ]
            print('\t'.join([str(o) for o in out]), file=open_output)

    pg_conn.close()

def main():
    argparse = ArgumentParser()
    argparse.add_argument('--input', help='Path to the file containing the taxonomies and assemblies', required=True)
    argparse.add_argument('--properties_dir', help='Path to the directory where the release1 application.properties are stored', required=True)
    argparse.add_argument('--assembly_dirs',  help='Path to the directory containing pre-downloaded species assemblies', required=True, nargs='+')
    argparse.add_argument('--download_dir', help='Path to the temporary directory where additional species assemblies will be downloaded', required=True,)
    argparse.add_argument('--release2_reference_folder', help='Path to the directory where selected fasta and report will be copied', required=True)
    argparse.add_argument('--output_assemblies_tsv', help='Path to the tsv file that will contain the list of assemblies to process', required=True)
    argparse.add_argument('--output_taxonomy_tsv', help='Path to the tsv file that will contain the list of species to process', required=True)
    argparse.add_argument('--eva_accession_path', help='path to the directory that contain eva-accession code and private json file.')
    argparse.add_argument("--private_config_xml_file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    argparse.add_argument('--debug', help='Set login level to debug', action='store_true', default=False)

    args = argparse.parse_args()
    log_cfg.add_stdout_handler()
    if args.debug:
        log_cfg.set_log_level(level=logging.DEBUG)

    global eva_accession_path
    if args.eva_accession_path:
            eva_accession_path = args.eva_accession_path
    
    aggregate_list_of_species(args.input, args.properties_dir, args.assembly_dirs, args.download_dir, args.release2_reference_folder,
                              args.output_assemblies_tsv, args.output_taxonomy_tsv, args.private_config_xml_file)


if __name__ == "__main__":
    main()
