#!/usr/bin/env python
import csv
from argparse import ArgumentParser
from collections import defaultdict


def parse_input(input_file):
    """
    Parse the input with a dict reader and returns the row organised per EVA/dbSNP then taxonomy id then assembly id
    """
    dbsnp_data = defaultdict(dict)
    eva_data = defaultdict(dict)
    
    with open(input_file) as open_file:
        reader = csv.DictReader(open_file, delimiter='\t')
        for row in reader:
            if row['Source'] == 'DBSNP':
                dbsnp_data[row['Taxid']][row['Assembly']] = row
            else:
                eva_data[row['Taxid']][row['Assembly']] = row
    return eva_data, dbsnp_data


def process_eva_assemblies(eva_data, dbsnp_data):
    """
    This method marks EVA assemblies as included or excluded. Will exclude them if they match any of the following rules.
     - A different assembly exists in dbSNP data for the same species (taxid)
     - The assembly used in EVA is from a different species.
     - Multiple EVA assemblies exists and this one is either not the one supported by Ensembl or has fewer variants.

     Only the first two rules currently excluded any assembly.
    """
    for taxid in eva_data:
        for assembly in eva_data[taxid]:
            row = eva_data.get(taxid).get(assembly)
            row['Decision to include'] = 'yes'
            row['Reason'] = ''
            if dbsnp_data.get(taxid, {}).get(assembly, {}) and len(dbsnp_data.get(taxid, {})) > 1 or \
                    assembly not in dbsnp_data.get(taxid, {}) and len(dbsnp_data.get(taxid, {})) > 0:
                # There is another assembly in dbSNP data that does not match the current one.
                # That could cause the variants from EVA to not be cluster with existing dbSNP RSid
                # which in turn could cause large number of variants to be retracted next release.
                row['Reason'] = 'Conflicting Assembly in dbSNP for this species'
                row['Decision to include'] = 'no'
            elif taxid != row.get('Taxid From Assembly'):
                # The assembly used to align these variants is not from the same taxonmy id.
                # If we cluster these variants it could "contaminate" this assembly's variants
                # Until we label them properly we cannot incorporate them.
                row['Reason'] = 'Species uses assembly from a different species.'
                row['Decision to include'] = 'no'

        # Count the number of assembly left
        assemblies_left = [
            local_assembly
            for local_assembly in eva_data.get(taxid)
            if eva_data[taxid][local_assembly]['Decision to include'] == 'yes'
        ]
        if len(assemblies_left) > 1:
            # We cannot release both assemblies so select the one associated with which supported by Ensembl.
            # If none are supported then release the one that has the most variants
            ensembl_associated_assembly = [
                local_assembly
                for local_assembly in assemblies_left
                if eva_data[taxid][local_assembly].get('Ensembl Assembly From Taxid') == local_assembly
            ]
            if ensembl_associated_assembly:
                assembly_to_keep = ensembl_associated_assembly[0]
                reason = 'Multiple EVA assemblies, can only choose one, so chose one that match Ensembl'
            else:
                assembly_to_keep = sorted(
                    assemblies_left,
                    key=lambda local_assembly: int(eva_data[taxid][local_assembly]['Number Of Variants']),
                )[-1]
                reason = 'Multiple EVA assemblies, can only choose one, so chose one with most variants'

            assemblies_left.remove(assembly_to_keep)
            for assembly in assemblies_left:
                row = eva_data.get(taxid).get(assembly)
                row['Reason'] = reason
                row['Decision to include'] = 'no'


def process_dbsnp_assemblies(eva_data, dbsnp_data):
    """
    This method marks dbSNP assemblies as included
    """
    for taxid in dbsnp_data:
        for assembly in dbsnp_data[taxid]:
            row = dbsnp_data.get(taxid).get(assembly)
            row['Decision to include'] = 'yes'
            row['Reason'] = 'Was released previously an need to be released again'


def write_output(output_file, eva_data, dbsnp_data):
    output_headers = [
        'Source', 'Assembly', 'Taxid', 'Number Of Studies', 'Number Of Variants', 
        'Taxid From Assembly', 'Scientific Name From Assembly', 'Ensembl Assembly From Assembly', 
        'Taxid From Taxid', 'Scientific Name From Taxid', 'Ensembl Assembly From Taxid',
        'Decision to include', 'Reason'
    ]
    with open(output_file, 'w') as open_file:
        writer = csv.DictWriter(open_file, fieldnames=output_headers, delimiter='\t')
        writer.writeheader()
        for taxid in eva_data:
            for assembly in eva_data[taxid]:
                writer.writerow(eva_data[taxid][assembly])
        for taxid in dbsnp_data:
            for assembly in dbsnp_data[taxid]:
                writer.writerow(dbsnp_data[taxid][assembly])


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--input', help='Path to the file containing the taxonomies and assemblies', required=True)
    argparse.add_argument('--output', help='Path to the file that will contain the input plus annotation', required=True)

    args = argparse.parse_args()
    
    eva_data, dbsnp_data = parse_input(args.input)
    process_eva_assemblies(eva_data, dbsnp_data )
    process_dbsnp_assemblies(eva_data, dbsnp_data)
    write_output(args.output, eva_data, dbsnp_data)


if __name__ == "__main__":
    main()
