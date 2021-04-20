import glob
import os
import tempfile
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from tasks.eva_2389.run_vcf_vertical_concat_pipeline import run_vcf_vertical_concat_pipeline, get_output_vcf_file_name
from unittest import TestCase


class TestVCFVerticalConcat(TestCase):
    # Tests require nextflow and bcftools installed locally and in PATH
    def test_concat(self):
        with tempfile.TemporaryDirectory() as tempdir:
            vcf_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "resources")
            #   s0.vcf.gz   s1.vcf.gz   s2.vcf.gz   s3.vcf.gz   s4.vcf.gz
            #       \           /           \           /
            #        s01.vcf.gz               s23.vcf.gz        s4.vcf.gz       <-------- Stage 0
            #               \                      /
            #                \                   /
            #                 \                /
            #                   s0123.vcf.gz                    s4.vcf.gz       <-------- Stage 1
            #                           \                       /
            #                            \                    /
            #                             \                 /
            #                                 Final_merged                      <-------- Stage 2
            run_vcf_vertical_concat_pipeline(toplevel_vcf_dir=vcf_dir, concat_processing_dir=tempdir,
                                             concat_chunk_size=2, bcftools_binary="bcftools",
                                             nextflow_binary="nextflow", nextflow_config_file=None, resume=False)
            stage_dirs = glob.glob(f"{tempdir}/vertical_concat/stage*")
            self.assertEqual(3, len(stage_dirs))
            output_vcf_from_multi_stage_concat = get_output_vcf_file_name(concat_stage_index=2, concat_batch_index=0,
                                                                          concat_processing_dir=tempdir)

            input_vcfs = sorted(glob.glob(f"{vcf_dir}/*.vcf.gz"))
            output_vcf_from_single_stage_concat = f"{tempdir}/single_stage_concat_result.vcf.gz"
            run_command_with_output("Concatenate VCFs with single stage...", f"bcftools concat {' '.join(input_vcfs)} "
                                                                             f"--allow-overlaps --remove-duplicates "
                                                                             f"-O z "
                                                                             f"-o {output_vcf_from_single_stage_concat}"
                                    )
            diffs = run_command_with_output("Compare outputs from single and multi-stage concat processes...",
                                            f'bash -c "diff '
                                            f'<(zcat {output_vcf_from_single_stage_concat} | grep -v ^#) '
                                            f'<(zcat {output_vcf_from_multi_stage_concat} | grep -v ^#)"',
                                            return_process_output=True)
            self.assertEqual("", diffs.strip())
