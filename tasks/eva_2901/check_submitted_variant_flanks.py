import argparse
import os
import itertools

from Bio import pairwise2

from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.config_utils import get_primary_mongo_creds_for_profile
from ebi_eva_common_pyutils.taxonomy.taxonomy import get_scientific_name_from_ensembl
from pymongo import MongoClient

cache = {'scientific_name_from_taxonomy': {}}


def get_genome(taxonomy, assembly_accession):
    if taxonomy not in cache['scientific_name_from_taxonomy']:
        species_scientific_name = get_scientific_name_from_ensembl(taxonomy)
        cache['scientific_name_from_taxonomy'][taxonomy] = os.path.join(
            cfg.query('genome_downloader', 'output_directory'),
            species_scientific_name.lower().replace(' ', '_'),
            assembly_accession, assembly_accession + '.fa'
        )
    return cache['scientific_name_from_taxonomy'][taxonomy]


def compare_variant_flanks(sequence1, sequence2):
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    return max((a.score for a in alignments))


def format_output(variant1, variant2, score):
    print(variant1)
    print(variant2)
    print(score)


def check_submitted_variant_flanks(mongo_client, ssid):
    samtools = cfg.query('executable', 'samtools', ret_default='samtools')
    sve_collection = mongo_client['eva_accession_sharded']['dbsnpSubmittedVariantEntity']
    cursor = sve_collection.find({'accession': ssid, 'remappedFrom': {'$exists': False}})
    flank_size = 50
    variant_records = list(cursor)
    id_2_info = {}
    for variant_rec in variant_records:
        flank_up_coord = f"{variant_rec['contig']}:{variant_rec['start'] - flank_size}-{variant_rec['start'] - 1}"
        flank_down_coord = f"{variant_rec['contig']}:{variant_rec['start'] + 1}-{variant_rec['start'] + flank_size}"

        position = f"{variant_rec['contig']}:{variant_rec['start']}"
        change = f"{variant_rec['ref']}-{variant_rec['alt']}"
        genome_assembly_fasta = get_genome(assembly_accession=variant_rec['seq'])
        commands = f"{samtools} faidx {genome_assembly_fasta} {flank_up_coord} | grep -v '^>' | sed 's/\n//' "
        flank_up = run_command_with_output(f'Extract upstream sequence using {flank_up_coord}',  commands).upper()
        commands = f"{samtools} faidx {genome_assembly_fasta} {flank_down_coord} | grep -v '^>' | sed 's/\n//' "
        flank_down = run_command_with_output(f'Extract upstream sequence using {flank_down_coord}',  commands).upper()
        id_2_info[variant_rec['_id']] = {'variant_rec': variant_rec, 'flank_up': flank_up, 'flank_down': flank_down}

    for variant_id1, variant_id2 in list(itertools.combinations(id_2_info, 2)):
        score = compare_variant_flanks(
            id_2_info[variant_id1]['flank_up'] + id_2_info[variant_id1]['variant_rec']['ref'] + id_2_info[variant_id1]['flank_down'],
            id_2_info[variant_id2]['flank_up'] + id_2_info[variant_id1]['variant_rec']['ref'] + id_2_info[variant_id2]['flank_down']
        )
        format_output(variant_id1, variant_id2, score)


def load_config(*args):
    cfg.load_config_file(
        os.path.expanduser('~/.submission_config.yml'),
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--profile_name",
                        help="Maven profile to use when connecting to mongodb", default='development')
    parser.add_argument("--ssid", help="single ssid to query",  required=True)
    args = parser.parse_args()
    load_config()
    mongo_host, mongo_user, mongo_pass = get_primary_mongo_creds_for_profile(args.profile_name)
    mongo_uri = f'mongodb://{mongo_user}:@{mongo_host}:27017/admin'
    mongo_client = MongoClient(mongo_uri, password=mongo_pass)
    check_submitted_variant_flanks(mongo_client, args.ssid)
    mongo_client.close()
