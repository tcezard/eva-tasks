#!/usr/bin/env python
import csv
import glob
import os
from argparse import ArgumentParser
from collections import defaultdict

from tasks.eva_2121.assembly_downloader import download_assembly

total_number_variant = 0


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


def find_reference_and_assembly_report(species_name, accession, assembly_dir):
    fasta, report = None, None
    accession_dir = os.path.join(
        assembly_dir,
        species_name.lower().replace(' ', '_'),
        accession
    )
    os.makedirs(accession_dir, exist_ok=True)
    fasta, report = download_assembly(assembly_dir,accession)
    # fastas = glob.glob(os.path.join(accession_dir, '*.fa'))
    # nucleotide_fastas = glob.glob(os.path.join(accession_dir, '*.fna'))
    # reports = glob.glob(os.path.join(accession_dir, '*_assembly_report.txt'))
    # if fastas and len(fastas) == 1:
    #     fasta = fastas[0]
    # elif nucleotide_fastas and len(nucleotide_fastas) == 1:
    #     fasta = nucleotide_fastas[0]
    # else:
    #     print('Found %s fasta file with pattern: %s' % (len(fastas) + len(nucleotide_fastas), os.path.join(accession_dir, '*.f?a')))
    # if not reports:
    #     report, fasta = download_assembly(detination_dir=accession_dir, assembly_accession=accession)
    # elif reports:
    #     report = reports[0]
    # if report:
    #     print('Found %s assembly report file with pattern: %s' % (len(reports), os.path.join(accession_dir, '*_assembly_report.txt')))
    return fasta, report


def assign_tempmongo_host(number_variants):
    global total_number_variant
    total_number_variant += number_variants
    instance_number = (total_number_variant // 180000000) + 1
    return 'tempmongo-' + str(instance_number)


def aggregate_list_of_species(input_file, assembly_dir):
    """
    """
    data_per_taxid, data_per_taxid_and_assembly = parse_input(input_file)
    only_unmapped_tax_id = []
    unchanged_taxid = []
    dbsnps_only = []
    eva_only = []
    mixture = []
    for taxid in data_per_taxid:
        rows = data_per_taxid[taxid]
        sources = {row['Source'] for row in rows}
        stay_unchanged = {row['Will stay unchanged'] for row in rows}
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

    print(len(only_unmapped_tax_id), only_unmapped_tax_id)
    print(len(unchanged_taxid), unchanged_taxid)
    print(len(dbsnps_only), dbsnps_only)
    print(len(eva_only), eva_only)
    print(len(mixture), mixture)

    for taxid, assembly in data_per_taxid_and_assembly:
        if taxid not in mixture + eva_only + dbsnps_only:
            continue
        rows = data_per_taxid_and_assembly[(taxid, assembly)]
        scientific_name = {row['Scientific Name From Taxid'] for row in rows}
        assert len(scientific_name) == 1
        scientific_name = scientific_name.pop()
        if assembly != 'Unmapped':
            fasta, report = find_reference_and_assembly_report(scientific_name, assembly, assembly_dir)
        if not report:
            fasta, report = None, None
        number_variants = sum([int(row['Number Of Variants (submitted variants)'].replace(',', '')) for row in rows])
        out = [
            taxid,
            scientific_name,
            assembly,
            ','.join({row['Source'] for row in rows}),
            fasta,
            report,
            assign_tempmongo_host(number_variants)
        ]
        print('\t'.join([str(o) for o in out]))


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--input', help='Path to the file containing the taxonomies and assemblies', required=True)
    argparse.add_argument('--assembly_dir', help='Path to the directory containing all species assemblies', required=True)
    args = argparse.parse_args()
    
    aggregate_list_of_species(args.input, args.assembly_dir)


if __name__ == "__main__":
    main()
