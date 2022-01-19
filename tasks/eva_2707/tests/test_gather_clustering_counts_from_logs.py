import os
from unittest import TestCase

from tasks.eva_2707.gather_clustering_counts_from_logs import parse_one_log


class TestGatherClusteringCounts(TestCase):
    steps = [
        'CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP', 'PROCESS_RS_MERGE_CANDIDATES_STEP',
        'PROCESS_RS_SPLIT_CANDIDATES_STEP', 'CLUSTERING_NON_CLUSTERED_VARIANTS_FROM_MONGO_STEP'
    ]

    metrics = [
        'clustered_variants_created',
        'clustered_variants_updated',
        'clustered_variants_rs_split',
        'clustered_variants_merged',
        'clustered_variants_merge_operations',
        'submitted_variants_clustered',
        'submitted_variants_kept_unclustered',
        'submitted_variants_updated_rs',
        'submitted_variants_update_operations',
        'item_written'
    ]

    def test_parse_one_log_simple(self):
        log_file = os.path.join(os.path.dirname(__file__), 'cluster_20211111225415.log')
        results = parse_one_log(log_file)
        assert set(results.keys()) == set(self.steps)
        for step in self.steps:
            assert set(results[step].keys()) == set(self.metrics)
            print(results[step].items())

    def test_parse_one_log_partial(self):
        # Log file for Maize where the second process didn't run
        steps = ['CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP']
        log_file = os.path.join(os.path.dirname(__file__), 'cluster_20211219230839.log')
        results = parse_one_log(log_file)
        assert set(results.keys()) == set(steps)
        for step in steps:
            assert set(results[step].keys()) == set(self.metrics)

    def test_parse_one_log_truncated(self):
        # Log file for horse where crash occured during the run which resulted in no count being logged
        step = 'CLUSTERING_CLUSTERED_VARIANTS_FROM_MONGO_STEP'
        log_file = os.path.join(os.path.dirname(__file__), 'cluster_20211125213810.log')
        results = parse_one_log(log_file)
        assert list(results.keys()) == [step]
        # Only item_written is present because the counts were not logged
        assert results[step] == {'item_written': '13463600'}

