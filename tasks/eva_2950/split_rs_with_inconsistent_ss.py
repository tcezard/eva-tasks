import os
import ast
import glob
import json
from argparse import ArgumentParser
from collections import defaultdict

from ebi_eva_common_pyutils.logger import logging_config
from itertools import combinations
import datetime

from pyfaidx import Fasta

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()
# Code from VCFtidy (https://github.com/quinlan-lab/vcftidy/blob/master/vcftidy.py)


def leftalign(chrom, pos, ref, alt, fa, max_shift=1000):
    # Add the context base if it is not there already
    if ref == '':
        ref = fa[chrom][pos - 2]  # -1 because pos is 1-based and -1 because we want the base before
        pos -= 1
        alt = ref + alt
    if alt == '':
        alt = fa[chrom][pos - 2]
        pos -= 1
        ref = alt + ref
    seq = fa[chrom][max(0, pos - max_shift - 1):pos + len(ref) - 1]
    assert seq.endswith(ref), (chrom, pos, ref, alt, seq[-10:])
    return _leftalign(pos, ref, alt, seq)[:3]


def _leftalign(pos, ref, alt, seq):
    """
    simple implementation from the vt paper:
    # actual variant is 2-base CA insertion.
    Last argument indicates whether we ran out of sequence and therefore did not
    finish left-aligning before running out of sequence. (False is bad).
    >>> _leftalign(123, 'CAC', 'C', 'GGGCACACAC')
    (118, 'GCA', 'G', True)
    # run out of sequence!
    >>> _leftalign(123, 'CAC', 'C', 'CACACAC')
    (119, 'CAC', 'C', False)
    >>> _leftalign(123, 'CCA', 'CAA', 'ACCCCCCA')
    (123, 'CC', 'CA', True)
    # have to left-trim after left-align
    >>> normalize(*_leftalign(123, 'CCA', 'CAA', 'ACCCCCCA')[:3], left_only=True)
    (124, 'C', 'A')
    >>> _leftalign(123, 'C', 'A', 'ACCCCCC')
    (123, 'C', 'A', True)
    """
    assert seq.endswith(ref)
    assert ref != alt
    seq = seq[:-len(ref)]
    ref, alt = list(ref), list(alt)
    j = 0

    quit = False
    while j < len(seq) and not quit:
        quit = True

        if ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            quit = False

        if len(ref) == 0 or len(alt) == 0:
            j += 1
            ref = [seq[-j]] + ref
            alt = [seq[-j]] + alt
            quit = False

    return pos - j, "".join(ref), "".join(alt), quit


def normalize(pos, ref, alt, left_only=False):
    """simplify a ref/alt a la vt normalize so that ref=CC, alt=CCT becomes
    ref=C, alt=CT. this helps in annotating variants.
    This code relies on the observation by Eric Minikel that many annotation
    misses can be addressed by removing common suffix and prefixes.
    (http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/)
    >>> normalize(123, 'T', 'C')
    (123, 'T', 'C')
    >>> normalize(123, 'CC', 'CCT')
    (124, 'C', 'CT')
    >>> normalize(123, 'TCCCCT', 'CCCCT')
    (123, 'TC', 'C')
    >>> normalize(123, 'TCCCCTA', 'CCCCT')
    (123, 'TCCCCTA', 'CCCCT')
    >>> normalize(123, 'TCCCCTA', 'CCCCTA')
    (123, 'TC', 'C')
    >>> normalize(123, 'AAATCCCCTA', 'AAACCCCTA')
    (125, 'AT', 'A')
    >>> normalize(123, 'CAAATCCCCTAG', 'AAACCCCTA')
    (123, 'CAAATCCCCTAG', 'AAACCCCTA')
    """
    if len(ref) == len(alt) == 1:
        return pos, ref, alt

    # logic for trimming from either end is the same so we just reverse the
    # string to trim the suffix (i == 0) and correct position when doing prefix
    # (i == 1). To support alleles that have already been right-trimmed from
    # _leftalign, we allow the left_only argument to just do prefix-trimming.
    if left_only:
        sides = (1,)
    else:
        sides = (0, 1)
    for i in sides:
        if i == 0: # suffix so flip
            ref, alt = ref[::-1], alt[::-1]

        n, min_len = 0, min(len(ref), len(alt))
        while n + 1 < min_len and alt[n] == ref[n]:
            n += 1

        alt, ref = alt[n:], ref[n:]
        if i == 0: # flip back
            ref, alt = ref[::-1], alt[::-1]
        else: # add to position since we stripped the prefix
            pos += n

    return pos, ref, alt


def leftnorm(chrom, pos, ref, alt, fa=None):
    """
    this is the normalization function that should be used.
    if no fasta is present, then it just normalizes. Otherwise
    it left-aligns and then normalizes.
    """
    if fa is None:
        return normalize(pos, ref, alt)

    return normalize(*leftalign(chrom, pos, ref, alt, fa), left_only=True)


def ast_parse_with_datetime(astr):
    """Based on https://stackoverflow.com/questions/4235606/way-to-use-ast-literal-eval-to-convert-string-into-a-datetime"""
    try:
        tree = ast.parse(astr)
    except SyntaxError:
        raise ValueError(astr)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.Expr, ast.Dict, ast.Str, ast.List, ast.Constant,
                            ast.Attribute, ast.Num, ast.Name, ast.Load, ast.Tuple)):
            continue
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'datetime'):
            continue
        raise ValueError(astr)
    return eval(astr)


def parse_eva2850_diagnositc_log(log_file):
    rsid = None
    with open(log_file) as open_file:
        for line in open_file:
            line = line.strip()
            if rsid:
                # We are parsing the line just after finding the RS id
                list_of_ssids = ast_parse_with_datetime(line.strip()[11:])
                yield rsid, list_of_ssids
                rsid = None
            if 'Not all original SS has same info' in line:
                sp_line = line.split()
                rsid = int(sp_line[4].strip(','))

        if rsid:
            list_of_ssids = ast_parse_with_datetime(line.strip()[11:])
            yield rsid, list_of_ssids


genome_cache = {}


def get_genome_object(ref_genome_directory, genome_accession):
    if genome_accession not in genome_cache:
        path_to_search = os.path.join(ref_genome_directory, '*', genome_accession, genome_accession + '.fa')
        genome_paths = glob.glob(path_to_search)
        if len(genome_paths) == 0:
            raise ValueError(f'Cannot locate genome for {genome_accession} in {ref_genome_directory}')
        genome_path = genome_paths[0]
        genome_cache[genome_accession] = Fasta(genome_path, as_raw=True, read_ahead=40000)
    return genome_cache[genome_accession]


def process_diagnostic_log(log_file, ref_genome_directory):
    for rsid, list_of_ss_entities in parse_eva2850_diagnositc_log(log_file):
        variant_to_entities = defaultdict(list)
        for ss_entity in list_of_ss_entities:
            genome_object = get_genome_object(ref_genome_directory, ss_entity['seq'])
            normalised_variant_definition = leftnorm(
                ss_entity['contig'], ss_entity['start'], ss_entity['ref'], ss_entity['alt'], fa=genome_object
            )
            ss_entity['contig'] + str(ss_entity['start']) + ss_entity['ref'] + ss_entity['alt']

            variant_to_entities[normalised_variant_definition].append(ss_entity)
        if len(variant_to_entities) == 1:
            logger.info(f'No split required for rs{rsid}')
        else:
            logger.info(f'rs{rsid} must be split in {len(variant_to_entities)} ')


def main():
    parser = ArgumentParser()
    parser.add_argument('--diagnostic_file')
    parser.add_argument('--ref_genome_directory')
    args = parser.parse_args()

    process_diagnostic_log(args.diagnostic_file, args.ref_genome_directory)


if __name__ == '__main__':
    main()
