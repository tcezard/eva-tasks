import sys

import builtins
import os
import hashlib
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from unittest import TestCase
from unittest.mock import patch

from tasks.eva_2778.check_rs_exist import find_rs_entity_not_exist_in_collection, find_rs_references_in_ss_collection, \
    check_rs_for_assembly


def calculate_id(rs):
    h = hashlib.sha1()
    h.update(str(rs).encode())
    return h.hexdigest().upper()


class TestCheckRsExist(TestCase):
    def setUp(self) -> None:
        host = 'localhost'
        port = 27017

        # Accessioning warehouse
        self.accession_db = 'eva_accession_sharded'
        self.submitted_variants_collection = 'submittedVariantEntity'
        self.clustered_variants_collection = 'clusteredVariantEntity'
        self.dbsnp_clustered_variants_collection = 'dbsnpClusteredVariantEntity'
        uri = f'mongodb://{host}:{port}'
        self.mongo_db = MongoDatabase(uri=uri, db_name="eva_accession_sharded")
        self.connection_handle = self.mongo_db.mongo_handle

        # Accessioning db
        submitted_variants = [{
            "_id": calculate_id(rs),
            "seq": "GCA_000181335.4",
            "tax": 1111,
            "study": "PRJEB30318",
            "contig": "CM000001.1",
            "start": 76166296,
            "ref": "C",
            "alt": "T",
            "accession": 5318166021,
            "rs": rs,
        } for rs in range(1000, 1011)]

        clustered_variants = [{
            "_id": calculate_id(rs),
            "asm": "GCA_000181335.4",
            "tax": 1111,
            "contig": "CM000001.1",
            "start": 76166296,
            "accession": rs,
        } for rs in range(1000, 1010)]

        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle[self.accession_db][self.submitted_variants_collection].insert_many(submitted_variants)
        self.connection_handle[self.accession_db][self.clustered_variants_collection].drop()
        self.connection_handle[self.accession_db][self.clustered_variants_collection].insert_many(clustered_variants)

    def tearDown(self) -> None:
        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle[self.accession_db][self.clustered_variants_collection].drop()
        self.connection_handle.close()

    def test_find_rs_entity_not_exist_in_collection(self):
        missing_rsids = find_rs_entity_not_exist_in_collection(self.mongo_db, self.clustered_variants_collection,
                                                               rsid_list=[1000], assembly_accession='GCA_000181335.4')
        assert missing_rsids == []
        missing_rsids = find_rs_entity_not_exist_in_collection(self.mongo_db, self.clustered_variants_collection,
                                                               rsid_list=[1010], assembly_accession='GCA_000181335.4')
        assert missing_rsids == [1010]

    def test_find_rs_references_in_ss_collection(self):
        rs_lists = find_rs_references_in_ss_collection(self.mongo_db, self.submitted_variants_collection,
                                                       'GCA_000181335.4', batch_size=5)
        assert list(rs_lists) == [
            [1000, 1001, 1002, 1003, 1004],
            [1005, 1006, 1007, 1008, 1009],
            [1010]
        ]

    def test_check_rs_for_assembly(self):
        with patch('builtins.print') as mock_print:
            check_rs_for_assembly(self.mongo_db, 'GCA_000181335.4', batch_size=5)
            mock_print.assert_any_call('Found a EVA submitted variant entity referencing rs 1010 but no clustered '
                                       'variant entity was found for it.')
            mock_print.assert_any_call('Processes 5 EVA submitted variants', file=sys.stderr)
            mock_print.assert_any_call('Processes 10 EVA submitted variants', file=sys.stderr)
            mock_print.assert_any_call('Processes 11 EVA submitted variants', file=sys.stderr)
