import os
from unittest import TestCase
from tasks.eva_2083.decision_release2 import process_eva_assemblies


class TestDecisionRelease2(TestCase):

    species_list = os.path.join(os.path.dirname(os.path.abspath(__file__)), "combined_from_mongoprod_count_study_annotated.tsv")

    def test_process_eva_assemblies(self):

        # No corresponding assembly in dbSNP and Assembly taxonomy matches data's taxonomy
        eva_data = {'1000': {'asm1': {'Taxid From Assembly': '1000'}}}
        dbsnp_data = {}
        process_eva_assemblies(eva_data, dbsnp_data)
        self.assertEqual(eva_data['1000']['asm1']['Decision to include'], 'yes')

        # Matching assembly in dbSNP and Assembly taxonomy matches data's taxonomy
        eva_data = {'1000': {'asm1': {'Taxid From Assembly': '1000'}}}
        dbsnp_data = {'1000': {'asm1': {}}}
        process_eva_assemblies(eva_data, dbsnp_data)
        self.assertEqual(eva_data['1000']['asm1']['Decision to include'], 'yes')

        # Conflicting assembly in dbSNP and Assembly taxonomy matches data's taxonomy
        eva_data = {'1000': {'asm1': {'Taxid From Assembly': '1000'}}}
        dbsnp_data = {'1000': {'asm2': {}}}
        process_eva_assemblies(eva_data, dbsnp_data)
        self.assertEqual(eva_data['1000']['asm1']['Decision to include'], 'no')
        self.assertEqual(eva_data['1000']['asm1']['Reason'], 'Conflicting Assembly in dbSNP for this species')

        # No corresponding assembly in dbSNP but Assembly taxonomy does not match data's taxonomy
        eva_data = {'1000': {'asm1': {'Taxid From Assembly': '1001'}}}
        dbsnp_data = {}
        process_eva_assemblies(eva_data, dbsnp_data)
        self.assertEqual(eva_data['1000']['asm1']['Decision to include'], 'no')
        self.assertEqual(eva_data['1000']['asm1']['Reason'], 'Species uses assembly from a different species.')

        # Multiple EVA assemblies, first one matches Ensembl supported assembly
        eva_data = {'1000': {
            'asm1': {'Taxid From Assembly': '1000', 'Ensembl Assembly From Taxid': 'asm1', 'Number Of Variants': '10'},
            'asm2': {'Taxid From Assembly': '1000', 'Ensembl Assembly From Taxid': 'asm1', 'Number Of Variants': '11'}
        }}
        dbsnp_data = {}
        process_eva_assemblies(eva_data, dbsnp_data)
        self.assertEqual(eva_data['1000']['asm1']['Decision to include'], 'yes')
        self.assertEqual(eva_data['1000']['asm2']['Decision to include'], 'no')
        self.assertEqual(eva_data['1000']['asm2']['Reason'], 'Multiple EVA assemblies, can only choose one, so chose one that match Ensembl')

        # Multiple EVA assemblies, second ones has more variants
        eva_data = {'1000': {
            'asm1': {'Taxid From Assembly': '1000', 'Number Of Variants': '10'},
            'asm2': {'Taxid From Assembly': '1000', 'Number Of Variants': '11'}
        }}
        dbsnp_data = {}
        process_eva_assemblies(eva_data, dbsnp_data)
        self.assertEqual(eva_data['1000']['asm1']['Decision to include'], 'no')
        self.assertEqual(eva_data['1000']['asm1']['Reason'],
                         'Multiple EVA assemblies, can only choose one, so chose one with most variants')
        self.assertEqual(eva_data['1000']['asm2']['Decision to include'], 'yes')
