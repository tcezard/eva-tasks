import os
from unittest import TestCase

from tasks.eva_2301.find_FP_variants import process_files, assessment_result_per_variants_from_vcf


class TestFindFPVariants(TestCase):


    def test_process_files(self):
        ressource_dir = os.path.join(os.path.dirname(__file__), 'ressource')
        assessment_vcf_file = os.path.join(ressource_dir, 'NoBED_hap.py_assessement.vcf.gz')
        realigned_vcf_file = os.path.join(ressource_dir, 'HG002_GRCh37to38_fake_genotypes.vcf.gz')
        bam_file = os.path.join(ressource_dir, 'reads_aligned.sorted.bam')
        process_files(bam_file, realigned_vcf_file, assessment_vcf_file)
