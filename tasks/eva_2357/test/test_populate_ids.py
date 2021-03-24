from mock import patch
from pymongo import MongoClient
from unittest import TestCase
from tasks.eva_2357.populate_ids import populate_ids, get_assembly_report_rows


class TestCorrectChr(TestCase):
    def setUp(self) -> None:
        self.host = 'localhost'
        self.port = 27017
        # Variant warehouse
        self.variant_warehouse_db = 'eva_fcatus_90'
        self.variant_collection = 'variants_2_0'

        # Variant warehouse
        self.accession_db = 'eva_accession_test'
        self.submitted_variants_collection = 'submittedVariantEntity'

        self.connection_handle = MongoClient(self.host)

        # Insert to variant warehouse
        variant = {
            "_id": "NC_018728.3_76166296_C_T",
            "chr": "NC_018728.3",
            "start": 76166296,
            "files": [
                {
                    "fid": "ERZ795310",
                    "sid": "PRJEB30318",
                },
                {
                    "fid": "ERZ111111",
                    "sid": "PRJEB11111",
                }
            ],
            "type": "SNV",
            "ref": "C",
            "alt": "T",
            "ids": ["ss1"]
        }

        self.connection_handle[self.variant_warehouse_db][self.variant_collection].drop()
        self.connection_handle[self.variant_warehouse_db][self.variant_collection].insert_many([variant])

        # Insert to accessioning db
        submitted_variant = {
            "_id": "56D2BAFE38814E88A8F5127333034A31C96501DD",
            "seq": "GCA_000181335.4",
            "tax": 9685,
            "study": "PRJEB30318",
            "contig": "CM001383.3",
            "start": 76166296,
            "ref": "C",
            "alt": "T",
            "accession": "5318166021",
            "rs": "1000",
            "version": 1,
            "createdDate": "2019-07-08T07:42:41.492Z"
        }

        submitted_variant2 = {
            "_id": "EEA0111F7EBE1895A60E00B4495B3C607D9B4A36",
            "seq": "GCA_000181335.4",
            "tax": 9685,
            "study": "PRJEB11111",
            "contig": "CM001383.3",
            "start": 76166296,
            "ref": "C",
            "alt": "T",
            "accession": "2000",
            "rs": "1000",
            "version": 1,
            "createdDate": "2020-07-08T07:42:41.492Z"
        }

        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle[self.accession_db][self.submitted_variants_collection].insert_many([submitted_variant, submitted_variant2])

    def tearDown(self) -> None:
        self.connection_handle[self.variant_warehouse_db][self.variant_collection].drop()
        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle.close()

    @patch('tasks.eva_2357.populate_ids.get_mongo_uri_for_eva_profile')
    def test_correct(self, mock_get_mongo_uri_for_eva_profile):
        mock_get_mongo_uri_for_eva_profile.return_value = 'mongodb://127.0.0.1:27017'
        populate_ids('../test/settings.xml', '../test/databases.txt', profile='localhost',
                     mongo_accession_db=self.accession_db)
        variant = (self.connection_handle[self.variant_warehouse_db][self.variant_collection].find_one(
            {'_id': 'NC_018728.3_76166296_C_T'}))
        self.assertListEqual(variant['ids'], ['ss1', 'ss5318166021', 'rs1000', 'ss2000'])
