from argparse import ArgumentParser

import pysam


def process_files(vcf_files_list):
    """Extract Basic info from all the VCF files provided to populate a tsv file"""
    with open(vcf_files_list) as open_list:
        for file_name in open_list:
            sample_name = file_name.split('.')[0]
            with pysam.VariantFile(file_name.strip(), 'r') as vcf_in:
                for vcf_record in vcf_in:
                    variant_type='SNP' if len(vcf_record.ref) == 1 and len(vcf_record.alts[0]) == 1 else 'INDEL'
                    print(
                        ':'.join((vcf_record.chrom, str(vcf_record.pos), vcf_record.ref, vcf_record.alts[0]))
                        + f"\t{vcf_record.info['AF']}\t{variant_type}\t{sample_name}"
                    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--vcf_files_list')
    args = parser.parse_args()

    process_files(args.vcf_files_list)


