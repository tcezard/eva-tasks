from argparse import ArgumentParser

import pysam


def process_files(vcf_files_list, output_tsv):
    """Extract Basic info from all the VCF files provided to populate a tsv file"""

    with open(vcf_files_list) as open_list, open(output_tsv, 'w') as open_output:
        for file_name in open_list:
            sample_name = file_name.split('.')[0]
            with pysam.VariantFile(file_name.strip(), 'r') as vcf_in:
                for vcf_record in vcf_in:
                    variant_type='SNP' if len(vcf_record.ref) == 1 and len(vcf_record.alts[0]) == 1 else 'INDEL'
                    print(
                        ':'.join((vcf_record.chrom, str(vcf_record.pos), vcf_record.ref, vcf_record.alts[0]))
                        + f"\t{vcf_record.info['AF']}\t{variant_type}\t{sample_name}", file=open_output
                    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--vcf_files_list')
    parser.add_argument('--output_tsv')
    args = parser.parse_args()

    process_files(args.vcf_files_list, args.output_tsv)



