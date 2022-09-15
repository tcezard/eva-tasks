import os
import ast
import glob
from argparse import ArgumentParser
from collections import defaultdict
from itertools import zip_longest

from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
import datetime

from pyfaidx import Fasta
from pymongo import MongoClient

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


# Code from VCFtidy (https://github.com/quinlan-lab/vcftidy/blob/master/vcftidy.py)
def leftalign(chrom, pos, ref, alt, fa, max_shift=1000):
    # Add the context base if it is not there already
    if ref == '':
        ref = fa[chrom][pos - 2].upper()  # -1 because pos is 1-based and -1 because we want the base before
        pos -= 1
        alt = ref + alt
    if alt == '':
        alt = fa[chrom][pos - 2].upper()
        pos -= 1
        ref = alt + ref
    seq = fa[chrom][max(0, pos - max_shift - 1):pos + len(ref) - 1].upper()
    ref = ref.upper()
    alt = alt.upper()
    if not seq.endswith(ref):
        raise ReferenceError('The reference bases in the variant are different from the reference sequence ' + str((chrom, pos, ref, alt, seq[-10:])))
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


def add_contex_base(fa, chrom, pos, ref, alt):
    # Add the context base if it is not there already
    if ref == '':
        ref = fa[chrom][pos - 2].upper()  # -1 because pos is 1-based and -1 because we want the base before
        pos -= 1
        alt = ref + alt
    if alt == '':
        alt = fa[chrom][pos - 2].upper()
        pos -= 1
        ref = alt + ref
    return pos, ref, alt


def ast_parse_with_datetime(astr):
    """
    Parse a string representing a python data structure into a python data structure
    Based on https://stackoverflow.com/questions/4235606/way-to-use-ast-literal-eval-to-convert-string-into-a-datetime
    """
    try:
        tree = ast.parse(astr)
    except SyntaxError:
        raise ValueError(astr)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.Expr, ast.Dict, ast.Str, ast.List, ast.Constant,
                            ast.Attribute, ast.Num, ast.Name, ast.Load, ast.Tuple)):
            continue
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'datetime':
            continue
        raise ValueError(astr)
    return eval(astr)


def parse_eva2850_diagnostic_log(log_file):
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


def variant_type(ref, alt):
    if len(ref) == len(alt) == 1:
        return 'SNP'
    if len(ref) > 1 and len(alt) > 1:
        return 'MNP'
    if len(ref) > 1:
        return 'DEL'
    return 'INS'


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def shelve_submitted_variant_entities(mongo_handle, submitted_variant_ids):
    output_collection = 'eva2950_dbsnpSubmittedVariantEntity'
    batch_size = 1000
    for batch_sve_ids in grouper(submitted_variant_ids, batch_size):
        query_filter = {'_id': {'$in': batch_sve_ids}}
        documents = [d for d in mongo_handle["eva_accession_sharded"]['dbsnpSubmittedVariantEntity'].find(query_filter)]
        mongo_handle["eva_accession_sharded"][output_collection].insert(documents)
        nb_sve = mongo_handle["eva_accession_sharded"][output_collection].estimated_document_count()
        assert nb_sve == len(batch_sve_ids), 'Not all variants were transfer to the output collection ' + output_collection
        response = mongo_handle['dbsnpSubmittedVariantEntity'].delete_many(query_filter)
        assert response.deleted_count == len(batch_sve_ids), 'Not all variants were deleted from dbsnpSubmittedVariantEntity'


def process_diagnostic_log(log_file, ref_genome_directory, mongo_handle=None):
    count_normalisation = count_splits = 0
    all_submitted_variant_ids = []
    for rsid, list_of_ss_entities in parse_eva2850_diagnostic_log(log_file):
        try:

            variant_to_entities = defaultdict(list)
            for ss_entity in list_of_ss_entities:
                all_submitted_variant_ids.append(ss_entity['_id'])
                genome_object = get_genome_object(ref_genome_directory, ss_entity['seq'])
                context_pos, context_ref, context_alt = add_contex_base(
                    genome_object, ss_entity['contig'], ss_entity['start'], ss_entity['ref'], ss_entity['alt']
                )
                normalised_pos, normalised_ref, normalised_alt = leftnorm(
                    ss_entity['contig'], context_pos, context_ref, context_alt, fa=genome_object
                )
                context_entity = {'start': context_pos, 'ref': context_ref, 'alt': context_alt}
                normalised_clustered_variant_definition = (normalised_pos, variant_type(normalised_ref, normalised_alt))
                ss_entity['contig'] + str(ss_entity['start']) + ss_entity['ref'] + ss_entity['alt']
                normalised_entity = {'start': normalised_pos, 'ref': normalised_ref, 'alt': normalised_alt}
                variant_to_entities[normalised_clustered_variant_definition].append((ss_entity, context_entity, normalised_entity))
            count_normalisation += process_renormalisation(variant_to_entities)
            count_splits += process_split(rsid, variant_to_entities)
        except ReferenceError as e:
            logger.error(f'Cannot Process rs{rsid} because one of the ssid cannot be normalised')
            logger.error(str(e))
    logger.info(f'{count_normalisation} submitted variants need to be renormalised')
    logger.info(f'{count_splits} clustered variants need to be created')

    if mongo_handle:
        shelve_submitted_variant_entities(mongo_handle, all_submitted_variant_ids)


def process_split(rsid, variant_to_entities):
    count_split = 0
    if len(variant_to_entities) == 1:
        logger.info(f'No split required for rs{rsid}')
    else:
        logger.info(f'rs{rsid} must be split in {len(variant_to_entities)}')
        # smallest SSids keeps the rs
        cluster_priority = sorted(
            variant_to_entities,
            key=lambda x: min([ss_entity['accession'] for ss_entity, _, _ in variant_to_entities[x]]),
        )
        for ss_entity, context_entity, normalised_entity in variant_to_entities[cluster_priority[0]]:
            logger.info(f'ss{ss_entity["accession"]} rs field remain the same.')

        for cluster in cluster_priority[1:]:
            for ss_entity, context_entity, normalised_entity in variant_to_entities[cluster]:
                logger.info(f'ss{ss_entity["accession"]} rs field needs to be updated')
            count_split += 1
    return count_split


def process_renormalisation(variant_to_entities):
    count_normalisation = 0
    for clustered_variant_entity in variant_to_entities:
        normalised_pos, vtype = clustered_variant_entity
        logger.info(f'  * {vtype} at {normalised_pos}')
        for original_variant, context_entity, normalised_variant in variant_to_entities[clustered_variant_entity]:
            if context_entity['start'] != normalised_variant['start']:
                logger.info(f'    ** ss{original_variant["accession"]} variant change from '
                            f'{context_entity["start"]}-{normalised_variant["ref"]}-{normalised_variant["alt"]} to '
                            f'{normalised_variant["start"]}-{context_entity["ref"]}-{context_entity["alt"]}')
                count_normalisation += 1
            else:
                logger.info(f'    ** ss{original_variant["accession"]} does not need to be normalised')
    return count_normalisation


def main():
    parser = ArgumentParser()
    parser.add_argument('--diagnostic_file')
    parser.add_argument('--ref_genome_directory')
    parser.add_argument('--settings_xml_file')
    parser.add_argument('--profile', default='development')
    args = parser.parse_args()
    if args.settings_xml_file:
        mongo_uri = get_mongo_uri_for_eva_profile(args.profile, args.settings_xml_file)
        with MongoClient(mongo_uri) as mongo_handle:
            process_diagnostic_log(args.diagnostic_file, args.ref_genome_directory, mongo_handle)
    else:
        process_diagnostic_log(args.diagnostic_file, args.ref_genome_directory)


if __name__ == '__main__':
    main()
