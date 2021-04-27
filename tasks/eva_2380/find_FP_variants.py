import os
import sys
from argparse import ArgumentParser
from collections import defaultdict

import pysam


def fetch_bases(fasta, contig, start, length):
    """
    Returns a subsection from a specified FASTA contig. The start coordinate is 1-based.
    """
    zero_base_start = start - 1
    end = zero_base_start + length
    new_ref = fasta.fetch(reference=contig, start=zero_base_start, end=end)
    return new_ref


def assessment_result_per_variants_from_vcf(vcf_file):
    # Two samples in this VCF first the truth set then the query set
    print('Read Assessment VCF: ' + vcf_file, file=sys.stderr)
    assessment_result_per_variants = {}
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        for vcf_record in vcf_in:
            try:
                if vcf_record.samples['QUERY']['BD'] in ['FP']:
                    key = (vcf_record.chrom, vcf_record.pos)
                    assessment_result_per_variants[key] = (vcf_record.samples['QUERY']['BD'], vcf_record.samples['QUERY']['BVT'])
            except KeyError:
                pass
    return assessment_result_per_variants


def get_variant_id(assessment_result_per_variants, vcf_file):
    print('Read realigned VCF: ' + vcf_file, file=sys.stderr)
    variant_id_to_assessment_result = {}
    variant_id_to_alignment_result = {}
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        for vcf_record in vcf_in:
            # Check around the location because normalisation might not give exact coordinates
            keys = [
                (vcf_record.chrom, vcf_record.pos),
                (vcf_record.chrom, vcf_record.pos + 1),
                (vcf_record.chrom, vcf_record.pos - 1),
                (vcf_record.chrom, vcf_record.pos + 2),
                (vcf_record.chrom, vcf_record.pos - 2)
            ]
            for key in keys:
                if key in assessment_result_per_variants:
                    variant_id_to_assessment_result[vcf_record.id] = assessment_result_per_variants[key]
                    variant_id_to_alignment_result[vcf_record.id] = vcf_record
                    # once found don't continue
                    break

    return variant_id_to_assessment_result, variant_id_to_alignment_result


def get_variant_from_vcf(variant_ids, vcf_file):
    print('Read standard VCF: ' + vcf_file, file=sys.stderr)
    variant_id_to_record = {}
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        for vcf_record in vcf_in:
            if vcf_record.id in variant_ids:
                variant_id_to_record[vcf_record.id] = vcf_record
    return variant_id_to_record


def output_alignment_records(variant_id_to_assessment_result, bam_file):
    base, ext = os.path.splitext(bam_file)
    output_bam = base + '_annotated' + ext
    output_csv = base + '_annotated.csv'
    with pysam.AlignmentFile(bam_file, "rb") as bam_in, \
            pysam.AlignmentFile(output_bam, "wb", template=bam_in) as bam_out, \
            open(output_csv, 'w') as open_csv:
        for sam_record in bam_in:
            sp_name = sam_record.query_name.split('|')
            rs_id = sp_name[4]
            if rs_id in variant_id_to_assessment_result:
                bd_tag, bvt_tag = variant_id_to_assessment_result[rs_id]
                sam_record.set_tag('BD', bd_tag)
                sam_record.set_tag('BT', bvt_tag)
            else:
                sam_record.set_tag('BD', 'Filtered')
                sam_record.set_tag('BT', 'Unknown')
            bam_out.write(sam_record)

            print('\t'.join(
                [sam_record.query_name, str(sam_record.mapping_quality)] +
                [
                    str(sam_record.get_tag(tag) if sam_record.has_tag(tag) else 'NaN')
                    for tag in ['AS', 'XS', 'XN', 'XM', 'XG', 'NM', 'MD', 'BD', 'BT']
                ]),
                file=open_csv)


def get_alignment_records(variant_ids, bam_file):
    print('Read alignment BAM: ' + bam_file, file=sys.stderr)
    variant_id_to_bam_records = defaultdict(list)
    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        for sam_record in bam_in:
            sp_name = sam_record.query_name.split('|')
            rs_id = sp_name[4]
            if rs_id in variant_ids:
                variant_id_to_bam_records[rs_id].append(sam_record)
    return variant_id_to_bam_records


def process_files(realigned_vcf_file, assessment_vcf_file, source_vcf, truth_vcf, bam_files, old_ref, new_ref):
    assessment_result_per_variants = assessment_result_per_variants_from_vcf(assessment_vcf_file)
    variant_id_to_assessment_result, variant_id_to_alignment_result = get_variant_id(assessment_result_per_variants, realigned_vcf_file)
    variant_id_to_source_record = get_variant_from_vcf(variant_id_to_assessment_result, source_vcf)
    variant_id_to_truth_record = get_variant_from_vcf(variant_id_to_assessment_result, truth_vcf)
    alignment_record_lists = []
    for bam_file in bam_files:
        alignment_record_lists.append(get_alignment_records(variant_id_to_assessment_result, bam_file))

    flank_len = 100
    old_fasta = pysam.FastaFile(old_ref)
    new_fasta = pysam.FastaFile(new_ref)

    for variant_id in variant_id_to_assessment_result:
        source_var = variant_id_to_source_record[variant_id]
        truth_var = variant_id_to_truth_record[variant_id]

        print('################################')
        print(variant_id)
        print('In source:' + str(source_var))
        print('In truth:' + str(truth_var))
        print('In realignment:' + str(variant_id_to_alignment_result[variant_id]))
        for alignment_records in alignment_record_lists:
            for alignment_record in alignment_records[variant_id]:
                print('Alignments:' + alignment_record.to_string())
            print('----')

        old_ref_left = fetch_bases(old_fasta, source_var.chrom, source_var.pos - flank_len, flank_len)
        old_ref = fetch_bases(old_fasta, source_var.chrom, source_var.pos, len(source_var.ref))
        old_ref_right = fetch_bases(old_fasta, source_var.chrom, source_var.pos + len(source_var.ref), flank_len)

        new_ref_left = fetch_bases(new_fasta, truth_var.chrom, truth_var.pos - flank_len, flank_len)
        new_ref = fetch_bases(new_fasta, truth_var.chrom, truth_var.pos, len(truth_var.ref))
        new_ref_right = fetch_bases(old_fasta, truth_var.chrom, truth_var.pos + len(truth_var.ref), flank_len)

        print(old_ref_left + "*" + old_ref + "*" + old_ref_right)
        print(new_ref_left + "*" + new_ref + "*" + new_ref_right)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--alignment_bams', nargs='+')
    parser.add_argument('--realigned_vcf')
    parser.add_argument('--assessment_vcf')
    parser.add_argument('--source_vcf')
    parser.add_argument('--truth_vcf')
    parser.add_argument('--old_ref')
    parser.add_argument('--new_ref')
    args = parser.parse_args()

    process_files(args.realigned_vcf, args.assessment_vcf, args.source_vcf, args.truth_vcf, args.alignment_bams,
                  args.old_ref, args.new_ref)


