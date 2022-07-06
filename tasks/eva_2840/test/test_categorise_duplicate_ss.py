import random
from pprint import pprint

from itertools import cycle

import sys

import builtins
import os
import hashlib
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from unittest import TestCase
from unittest.mock import patch

from tasks.eva_2778.check_rs_exist import find_rs_entity_not_exist_in_collection, find_rs_references_in_ss_collection, \
    check_rs_for_assembly
from tasks.eva_2840.categorise_duplicate_ss import categorise_batch_duplicate_ss


def calculate_id(ss):
    h = hashlib.sha1()
    h.update(str(ss).encode())
    return h.hexdigest().upper()


class TestCategoriseSsExist(TestCase):
    def setUp(self) -> None:
        host = 'localhost'
        port = 27017

        # Accessioning warehouse
        self.accession_db = 'eva_accession_sharded'
        self.submitted_variants_collection = 'dbsnpSubmittedVariantEntity'
        uri = f'mongodb://{host}:{port}'
        self.mongo_db = MongoDatabase(uri=uri, db_name="eva_accession_sharded")
        self.connection_handle = self.mongo_db.mongo_handle

        # Accessioning db
        submitted_variants = [
            {
                "seq": "GCA_000181335.4",
                "tax": 1111,
                "study": "PRJEB30318",
                "contig": "CM000001.1",
                "start": random.randint(1, 1000000),
                "ref": "C",
                "alt": "T",
                "accession": ss,
                "rs": 5318166021,
            } for ss in [1000, 1000, 1001, 1001]
        ]
        for submitted_variant in submitted_variants:
            submitted_variant['_id'] = calculate_id(submitted_variant['start'])

        for i, sub_var in enumerate(submitted_variants):
            if i % 4 == 0:
                sub_var['allelesMatch'] = True
            if i % 4 == 1:
                sub_var['mapWeight'] = 3
            if i % 4 == 3:
                sub_var['remappedFrom'] = 'GCA_000181335.3'
        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle[self.accession_db][self.submitted_variants_collection].insert_many(submitted_variants)

    def tearDown(self) -> None:
        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle.close()

    def test_categorise_batch_duplicate_ss(self):
        c = self.connection_handle[self.accession_db][self.submitted_variants_collection].find()
        ssids = [1000, 1001]
        annotated_ssids = categorise_batch_duplicate_ss(self.mongo_db, ssids, assembly_accession='GCA_000181335.4')
        pprint(annotated_ssids)

