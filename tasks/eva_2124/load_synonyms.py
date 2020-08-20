import os
import wget
from urllib.parse import urlparse
from tasks.eva_2124.get_assembly_report_url import get_assembly_report_url


def load_synonyms_for_assembly(assembly_accession):
    """
    Reads an assembly report and loads dictionaries to map names to genbank.
    Returns 5 dictionaries: by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc
    Example usage:
    genbank = by_name['1']['genbank']
    """
    print('Searching assembly report for {} ...'.format(assembly_accession))
    url = get_assembly_report_url(assembly_accession)
    asm_report_file = os.path.basename(urlparse(url).path)
    downloaded_already = os.path.isfile(asm_report_file)

    if downloaded_already:
        print('Assembly report was already downloaded, skipping download ...')
    else:
        print('Downloading assembly report ...')
        download_file_from_ftp(url)

    print('Parsing assembly report ...')
    by_name = dict()
    by_assigned_molecule = dict()
    by_genbank = dict()
    by_refseq = dict()
    by_ucsc = dict()
    with open(asm_report_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            is_chromosome = columns[1] == 'assembled-molecule' and columns[3] == 'Chromosome'
            synonyms = {'name': columns[0],
                        'assigned_molecule': columns[2] if is_chromosome else None,
                        'is_chromosome': is_chromosome,
                        'genbank': columns[4],
                        'refseq': columns[6],
                        'is_genbank_refseq_identical': (columns[5] == '='),
                        'ucsc': columns[9] if columns[9] != 'na' else None}

            by_name[columns[0]] = synonyms
            by_genbank[columns[4]] = synonyms
            by_refseq[columns[6]] = synonyms

            if synonyms['assigned_molecule'] is not None:
                by_assigned_molecule[columns[2]] = synonyms

            if synonyms['ucsc'] is not None:
                by_ucsc[columns[9]] = synonyms

    print('Loaded chromosome synonyms for assembly {}'.format(assembly_accession))
    return by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc


def download_file_from_ftp(url):
    return wget.download(url)


