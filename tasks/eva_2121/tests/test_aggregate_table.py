from unittest import TestCase

from tasks.eva_2121 import aggregate_table


class TestAggregateTable(TestCase):

    def test_count_variants(self):
        rows = [
            {'Will stay unchanged': 'yes', 'Number Of Variants (submitted variants)': '1,000', 'Source': 'DBSNP'},
            {'Will stay unchanged': 'yes', 'Number Of Variants (submitted variants)': '1', 'Source': 'DBSNP - filesystem'},
            {'Will stay unchanged': 'no', 'Number Of Variants (submitted variants)': '2', 'Source': 'EVA'}
        ]
        assert aggregate_table.count_variants(rows, only_variant_to_process=True) == 1002
        assert aggregate_table.count_variants(rows, only_variant_to_process=False) == 1003

        rows = [
            {'Will stay unchanged': 'yes', 'Number Of Variants (submitted variants)': '1,000', 'Source': 'DBSNP'},
            {'Will stay unchanged': 'yes', 'Number Of Variants (submitted variants)': '1','Source': 'DBSNP - filesystem'},
        ]
        assert aggregate_table.count_variants(rows, only_variant_to_process=True) == 0
        assert aggregate_table.count_variants(rows, only_variant_to_process=False) == 1001

    def test_assign_tempmongo_host_round_robin(self):
        aggregate_table.number_of_tempmongo_instances = 4
        total_number_variants = 100
        assert aggregate_table.assign_tempmongo_host_round_robin(10, total_number_variants) == 'tempmongo-1'
        assert aggregate_table.assign_tempmongo_host_round_robin(14, total_number_variants) == 'tempmongo-1'
        assert aggregate_table.assign_tempmongo_host_round_robin(1, total_number_variants) == 'tempmongo-1'
        assert aggregate_table.assign_tempmongo_host_round_robin(10, total_number_variants) == 'tempmongo-2'
        assert aggregate_table.assign_tempmongo_host_round_robin(10, total_number_variants) == 'tempmongo-2'
        assert aggregate_table.assign_tempmongo_host_round_robin(15, total_number_variants) == 'tempmongo-3'
        assert aggregate_table.assign_tempmongo_host_round_robin(10, total_number_variants) == 'tempmongo-3'
        assert aggregate_table.assign_tempmongo_host_round_robin(10, total_number_variants) == 'tempmongo-4'
        assert aggregate_table.assign_tempmongo_host_round_robin(10, total_number_variants) == 'tempmongo-4'
        assert aggregate_table.assign_tempmongo_host_round_robin(10, total_number_variants) == 'tempmongo-4'
