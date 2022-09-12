import os
from unittest import TestCase

from tasks.eva_2933.collect_variant_types_from_gvf_files import collect_structural_variant_types


class TestCollectVariantTypesFromGVFFiles(TestCase):

    def test_collect_structural_variant_types(self):
        gvf_root_path = os.path.join(os.path.dirname(__file__), "data")
        variant_types = collect_structural_variant_types(gvf_root_path)
        self.assertSetEqual({'deletion', 'copy_number_variation', 'complex_structural_alteration', 'copy_number_loss',
                             'duplication', 'copy_number_gain'}, variant_types)
