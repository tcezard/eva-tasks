import os
from argparse import ArgumentParser

import pysam


def assessment_result_per_variants_from_vcf(vcf_file):
    # Two samples in this VCF first the truth set then the query set
    assessment_result_per_variants = {}
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        for vcf_record in vcf_in:
            try:
                if vcf_record.samples['QUERY']['BD'] in ['TP', 'FP']:
                    assessment_result_per_variants[(vcf_record.chrom, vcf_record.pos)] = (vcf_record.samples['QUERY']['BD'], vcf_record.samples['QUERY']['BVT'])
            except KeyError:
                pass
    return assessment_result_per_variants


def get_variant_id(assessment_result_per_variants, vcf_file):
    variant_id_to_assessment_result = {}
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        for vcf_record in vcf_in:
            if (vcf_record.chrom, vcf_record.pos) in assessment_result_per_variants:
                variant_id_to_assessment_result[vcf_record.id] = assessment_result_per_variants[(vcf_record.chrom, vcf_record.pos)]
    return variant_id_to_assessment_result


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


def process_files(bam_file, realigned_vcf_file, assessment_vcf_file):
    assessment_result_per_variants = assessment_result_per_variants_from_vcf(assessment_vcf_file)
    variant_id_to_assessment_result = get_variant_id(assessment_result_per_variants, realigned_vcf_file)
    output_alignment_records(variant_id_to_assessment_result, bam_file)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--alignment_bam')
    parser.add_argument('--realigned_vcf')
    parser.add_argument('--assessment_vcf')
    args = parser.parse_args()

    process_files(args.alignment_bam, args.realigned_vcf, args.assessment_vcf)


