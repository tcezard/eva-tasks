import argparse
import os
import sys

import itertools

from Bio import pairwise2
from Bio.Seq import Seq
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.config_utils import get_primary_mongo_creds_for_profile
from ebi_eva_common_pyutils.taxonomy.taxonomy import get_scientific_name_from_ensembl
from pymongo import MongoClient

cache = {'scientific_name_from_taxonomy': {}}


def revcomp(seq):
    return str(Seq(seq).reverse_complement())


def get_scientific_name_from_taxonomy(taxonomy):
    if taxonomy not in cache['scientific_name_from_taxonomy']:
        species_scientific_name = get_scientific_name_from_ensembl(taxonomy)
        cache['scientific_name_from_taxonomy'][taxonomy] =species_scientific_name
    return cache['scientific_name_from_taxonomy'][taxonomy]


def get_genome(taxonomy, assembly_accession):
    return os.path.join(
        cfg.query('genome_downloader', 'output_directory'),
        get_scientific_name_from_taxonomy(taxonomy).lower().replace(' ', '_'),
        assembly_accession, assembly_accession + '.fa'
    )


def compare_variant_flanks(sequence1, sequence2):
    alignments1 = pairwise2.align.globalms(sequence1, sequence2, 1, -3, -10, -10, one_alignment_only=True)
    alignments2 = pairwise2.align.globalms(sequence1, revcomp(sequence2), 1, -3, -10, -10, one_alignment_only=True)

    if alignments1 and alignments2:
        if alignments2[0].score > alignments1[0].score:
            return alignments2[0], '-'
        else:
            return alignments1[0], '+'
    elif alignments1:
        return alignments1, '+'
    elif alignments2:
        return alignments2, '-'


def format_output(ssid, variant1, variant2, alignment, strand, flank_up1, flank_down1, flank_up2, flank_down2):
    ref1 = ' '.join((flank_up1, variant1['ref'], flank_down1))
    ref_bases = variant2['ref'] if strand == '+' else revcomp(variant2['ref'])
    ref2 = ' '.join((flank_up2, ref_bases, flank_down2))
    alt1 = ' '.join((flank_up1, variant1['alt'], flank_down1))
    alt_bases = variant2['alt'] if strand == '+' else revcomp(variant2['alt'])
    alt2 = ' '.join((flank_up2, alt_bases, flank_down2))

    out = [ssid,
        variant1['seq'], variant1['contig'], variant1['start'], variant1['ref'], variant1['alt'],
        variant2['seq'], variant2['contig'], variant2['start'], variant2['ref'], variant2['alt'],
        alignment.score, strand, alignment.seqA, alignment.seqB, ref1, ref2, alt1, alt2
    ]
    return '\t'.join([str(s) for s in out])


def check_submitted_variant_flanks(mongo_client, ssid):
    samtools = cfg.query('executable', 'samtools', ret_default='samtools')
    sve_collection = mongo_client['eva_accession_sharded']['dbsnpSubmittedVariantEntity']
    cursor = sve_collection.find({'accession': int(ssid), 'remappedFrom': {'$exists': False}})
    flank_size = 50
    variant_records = list(cursor)
    id_2_info = {}
    for variant_rec in variant_records:
        flank_up_coord = f"{variant_rec['contig']}:{variant_rec['start'] - flank_size}-{variant_rec['start'] - 1}"
        flank_down_coord = f"{variant_rec['contig']}:{variant_rec['start'] + 1}-{variant_rec['start'] + flank_size}"
        genome_assembly_fasta = get_genome(assembly_accession=variant_rec['seq'], taxonomy=variant_rec['tax'])
        command = f"{samtools} faidx {genome_assembly_fasta} {flank_up_coord} | grep -v '^>' | sed 's/\\n//' "
        flank_up = run_command_with_output(f'Extract upstream sequence using {flank_up_coord}',  command, return_process_output=True).strip().upper()
        command = f"{samtools} faidx {genome_assembly_fasta} {flank_down_coord} | grep -v '^>' | sed 's/\\n//' "
        flank_down = run_command_with_output(f'Extract upstream sequence using {flank_down_coord}',  command, return_process_output=True).strip().upper()
        id_2_info[variant_rec['_id']] = {'variant_rec': variant_rec, 'flank_up': flank_up, 'flank_down': flank_down}

    for variant_id1, variant_id2 in list(itertools.combinations(id_2_info, 2)):
        alignment, strand = compare_variant_flanks(
            id_2_info[variant_id1]['flank_up'] + id_2_info[variant_id1]['variant_rec']['ref'] + id_2_info[variant_id1]['flank_down'],
            id_2_info[variant_id2]['flank_up'] + id_2_info[variant_id1]['variant_rec']['ref'] + id_2_info[variant_id2]['flank_down']
        )
        output = format_output(
            ssid, id_2_info[variant_id1]['variant_rec'], id_2_info[variant_id2]['variant_rec'], alignment, strand,
            id_2_info[variant_id1]['flank_up'], id_2_info[variant_id1]['flank_down'],
            id_2_info[variant_id2]['flank_up'], id_2_info[variant_id2]['flank_down']
        )
        print(output)


def load_config(*args):
    cfg.load_config_file(
        os.path.expanduser('~/.submission_config.yml'),
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--profile_name",
                        help="Maven profile to use when connecting to mongodb")
    parser.add_argument("--ssid_file", help="file containing a single ssid per line", type=str, required=True)
    args = parser.parse_args()
    # Get the config file loaded
    load_config()
    # Connect to mongodb database
    profile_name = args.profile_name or cfg['maven']['environment']
    mongo_host, mongo_user, mongo_pass = get_primary_mongo_creds_for_profile(profile_name, cfg['maven']['settings_file'])
    mongo_uri = f'mongodb://{mongo_user}:@{mongo_host}:27017/eva_accession_sharded?authSource=admin'
    mongo_client = MongoClient(mongo_uri, password=mongo_pass)
    nb_ssids = 0
    with open(args.ssid_file) as open_file:
        for line in open_file:
            nb_ssids += 1
            if nb_ssids % 1000 == 0 :
                print(f'Processed {nb_ssids} ssids', file=sys.stderr)
            check_submitted_variant_flanks(mongo_client, line.strip())
    mongo_client.close()
