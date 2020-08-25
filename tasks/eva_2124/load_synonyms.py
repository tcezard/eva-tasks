import os
import wget
from urllib.parse import urlparse
from tasks.eva_2124.get_assembly_report_url import get_assembly_report_url
import logging


def load_synonyms_for_assembly(assembly_accession, assembly_report_file=None):
    """
    Reads an assembly report and loads dictionaries to map names to genbank.
    The assembly report will be automatically downloaded if the parameter assembly_report_file is None.
    Returns 5 dictionaries: by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc
    Example usage:
    genbank = by_name['1']['genbank']
    """
    if assembly_report_file is None:
        assembly_report_file = download_assembly_report(assembly_accession)
    else:
        logging.info('Using provided assembly report at {}'.format(assembly_report_file))

    logging.info('Parsing assembly report ...')
    by_name = dict()
    by_assigned_molecule = dict()
    by_genbank = dict()
    by_refseq = dict()
    by_ucsc = dict()
    with open(assembly_report_file, 'r') as f:
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

    logging.info('Loaded chromosome synonyms for assembly {}'.format(assembly_accession))
    return by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc


def download_assembly_report(assembly_accession):
    logging.info('Searching assembly report for {} ...'.format(assembly_accession))
    url = get_assembly_report_url(assembly_accession)
    asm_report_file = os.path.basename(urlparse(url).path)
    downloaded_already = os.path.isfile(asm_report_file)
    if downloaded_already:
        logging.info('Assembly report was already downloaded, skipping download ...')
    else:
        logging.info('Downloading assembly report ...')
        download_file_from_ftp(url)
    return asm_report_file


def download_file_from_ftp(url):
    return wget.download(url)


