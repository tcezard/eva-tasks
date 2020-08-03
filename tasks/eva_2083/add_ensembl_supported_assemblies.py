#!/usr/bin/env python
import re
from argparse import ArgumentParser
from collections import defaultdict

import requests

eutils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
esearch_url = eutils_url + 'esearch.fcgi'
esummary_url = eutils_url + 'esummary.fcgi'
efetch_url = eutils_url + 'efetch.fcgi'
ensembl_url = 'http://rest.ensembl.org/info/assembly'

cache = defaultdict(dict)


def retrieve_species_names_from_tax_id(taxid):
    if taxid not in cache['taxid_to_name']:
        payload = {'db': 'Taxonomy', 'id': taxid}
        r = requests.get(efetch_url, params=payload)
        match = re.search('<Rank>(.+?)</Rank>', r.text, re.MULTILINE)
        rank = None
        if match:
            rank = match.group(1)
        if rank in ['species', 'subspecies']:
            scientific_name = None
            match = re.search('<ScientificName>(.+?)</ScientificName>', r.text, re.MULTILINE)
            if match:
                scientific_name = match.group(1)
                cache['taxid_to_name'][taxid] = scientific_name
        else:
            print('WARNING: No species found for %s' % taxid)
    return taxid, cache['taxid_to_name'].get(taxid)


def retrieve_species_name_from_assembly_accession(assembly_accession):
    if assembly_accession not in cache['assembly_to_species']:
        payload = {'db': 'Assembly', 'term': '"{}"'.format(assembly_accession), 'retmode': 'JSON'}
        data = requests.get(esearch_url, params=payload).json()
        if data:
            assembly_id_list = data.get('esearchresult').get('idlist')
            payload = {'db': 'Assembly', 'id': ','.join(assembly_id_list), 'retmode': 'JSON'}
            summary_list = requests.get(esummary_url, params=payload).json()
            all_species_names = set()
            for assembly_id in summary_list.get('result', {}).get('uids', []):
                assembly_info = summary_list.get('result').get(assembly_id)
                all_species_names.add((assembly_info.get('speciestaxid'), assembly_info.get('speciesname')))
            if len(all_species_names) == 1:
                cache['assembly_to_species'][assembly_accession] = all_species_names.pop()
            else:
                print('WARNING: %s taxons found for assembly %s ' % (len(all_species_names), assembly_accession))
    return cache['assembly_to_species'].get(assembly_accession) or (None, None)


def retrieve_current_ensembl_assemblies(source):
    """
    Retrieve the assembly accession currently supported by ensembl for the provided taxid or assembly accession
    In both case it looks up the associated species name in NCBI and using the species name returns the currently
    supported assembly for this species.
    """
    print('Search for species name for ' + source)
    if str(source).isdigit():
        taxid, scientific_name = retrieve_species_names_from_tax_id(source)
    else:
        taxid, scientific_name = retrieve_species_name_from_assembly_accession(source)

    if scientific_name:
        print('Found ' + scientific_name)
        if scientific_name not in cache['scientific_name_to_ensembl']:
            url = ensembl_url + '/' + scientific_name.lower().replace(' ', '_')
            response = requests.get(url, params={'content-type': 'application/json'})
            data = response.json()
            assembly_accession = str(data.get('assembly_accession'))
            cache['scientific_name_to_ensembl'][scientific_name] = assembly_accession
        return [str(taxid), str(scientific_name), cache['scientific_name_to_ensembl'].get(scientific_name)]

    return ['NA', 'NA', 'NA']


def main():
    argparse = ArgumentParser()
    argparse.add_argument('--input', help='Path to the file containing the taxonomies and assemblies', required=True)
    argparse.add_argument('--output', help='Path to the file that will contain the input plus annotation', required=True)

    args = argparse.parse_args()
    with open(args.input) as open_input, open(args.output, 'w') as open_ouput:
        for line in open_input:
            # Assume format of the input file is <assembly> <taxid> <...>
            sp_line = line.split()
            sp_line += retrieve_current_ensembl_assemblies(sp_line[0])
            sp_line += retrieve_current_ensembl_assemblies(sp_line[1])
            open_ouput.write('\t'.join(sp_line) + '\n')


if __name__ == "__main__":
    main()
