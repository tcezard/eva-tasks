import os.path
from pprint import pprint
from unittest import TestCase

from pyfaidx import Fasta

from tasks.eva_2950.split_rs_with_inconsistent_ss import leftnorm, parse_eva2850_diagnositc_log, process_diagnostic_log


class TestNormalisation(TestCase):

    test_dir = os.path.dirname(__file__)

    def test_leftnorm(self):
        fasta_file = os.path.join(self.test_dir, 'fasta_file.fa')
        # Based on rs54131737
        fa = Fasta(fasta_file, as_raw=True, read_ahead=40000)
        chrom = 'AP014957.1'
        pos = 91
        ref = ''
        alt = 'TTTTT'
        assert (leftnorm(chrom, pos, ref, alt, fa=fa)) == (80, 'C', 'CTTTTT')
        pos = 81
        ref = ''
        alt = 'TTTTT'
        assert (leftnorm(chrom, pos, ref, alt, fa=fa)) == (80, 'C', 'CTTTTT')
        pos = 81
        ref = ''
        alt = 'T'
        assert (leftnorm(chrom, pos, ref, alt, fa=fa)) == (80, 'C', 'CT')
        pos = 91
        ref = ''
        alt = 'T'
        assert (leftnorm(chrom, pos, ref, alt, fa=fa)) == (80, 'C', 'CT')


class TestSplitRS(TestCase):

    test_dir = os.path.dirname(__file__)

    def test_parse_eva2850_diagnositc_log(self):
        diagnostic_file = os.path.join(self.test_dir, 'diagnostic_output_log.out')
        rsids = []
        lists_of_ssids = []
        for rsid, list_of_ssids in parse_eva2850_diagnositc_log(diagnostic_file):
            rsids.append(rsid)
            lists_of_ssids.append(list_of_ssids)

        assert [l['accession'] for list_of_ssids in lists_of_ssids for l in list_of_ssids] == [
            71656146, 1961656906, 73429405, 71656146, 73429405, 73630724, 73526101, 1965862538, 73630724, 73526101,
            1964358410, 73565740, 1964358413, 73565740
        ]
        assert rsids == [54131737, 53378121, 54319631]

    def test_process_diagnostic_log(self):
        diagnostic_file = os.path.join(self.test_dir, 'diagnostic_output_log.out')
        ref_genome_dir = os.path.join(self.test_dir, 'references')

        process_diagnostic_log(diagnostic_file, ref_genome_dir)
