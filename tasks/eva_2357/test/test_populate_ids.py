import os

import pytest
from pymongo import MongoClient
from unittest import TestCase
from unittest.mock import patch

import logging
from tasks.eva_2357.populate_ids import populate_ids, check_all_contigs


class TestCorrectChr(TestCase):
    def setUp(self) -> None:
        self.host = 'localhost'
        self.port = 27017
        # Variant warehouse
        self.variant_warehouse_db = 'eva_fcatus_90'
        self.variant_warehouse_db_errors = 'eva_fcatus_90_errors'
        self.variant_collection = 'variants_2_0'

        # Variant warehouse
        self.accession_db = 'eva_accession_test'
        self.submitted_variants_collection = 'submittedVariantEntity'

        self.connection_handle = MongoClient(self.host)

        # Variant warehouse
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
                },
                {
                    "fid": "ERZ_NO_VARIANT",
                    "sid": "PRJEB_NO_VARIANT",
                }
            ],
            "type": "SNV",
            "ref": "C",
            "alt": "T",
            "ids": ["ss1"]
        }

        variant_not_in_acc_db = {
            "_id": "NC_018728.3_1000_G_T",
            "chr": "NC_018728.3",
            "start": 1000,
            "files": [
                {
                    "fid": "ERZ_0",
                    "sid": "PRJEB_0",
                }
            ],
            "type": "SNV",
            "ref": "G",
            "alt": "T"
        }

        variant_no_contig_synonym = {
            "_id": "NC_000001.1_2000_A_T",
            "chr": "NC_000001.1",
            "start": 2000,
            "files": [
                {
                    "fid": "ERZ_0",
                    "sid": "PRJEB_0",
                }
            ],
            "type": "SNV",
            "ref": "A",
            "alt": "T"
        }

        self.connection_handle[self.variant_warehouse_db][self.variant_collection].drop()
        self.connection_handle[self.variant_warehouse_db][self.variant_collection].insert_many([variant,
                                                                                                variant_not_in_acc_db])

        self.connection_handle[self.variant_warehouse_db_errors][self.variant_collection].drop()
        self.connection_handle[self.variant_warehouse_db_errors][self.variant_collection].insert_many([variant,
                                                                                                       variant_not_in_acc_db,
                                                                                                       variant_no_contig_synonym])

        # Accessioning db
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
        self.connection_handle[self.accession_db][self.submitted_variants_collection].insert_many([submitted_variant,
                                                                                                   submitted_variant2])

    def tearDown(self) -> None:
        self.connection_handle[self.variant_warehouse_db][self.variant_collection].drop()
        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle.close()

    @patch('tasks.eva_2357.populate_ids.get_mongo_uri_for_eva_profile')
    def test_populate_ids(self, mock_get_mongo_uri_for_eva_profile):
        logging.getLogger().setLevel(logging.DEBUG)
        mock_get_mongo_uri_for_eva_profile.return_value = 'mongodb://127.0.0.1:27017'
        settings = self.get_test_resource("settings.xml")
        databases_file = self.get_test_resource("databases.txt")
        self.assertEqual(1, populate_ids(settings, databases_file,  profile='localhost',
                                         mongo_accession_db=self.accession_db))
        variant = (self.connection_handle[self.variant_warehouse_db][self.variant_collection].find_one(
            {'_id': 'NC_018728.3_76166296_C_T'}))
        # The name for assertCountEqual more than just counting the elements, it checks both lists have the same
        # elements regardless of the order
        self.assertCountEqual(variant['ids'], ['ss1', 'ss5318166021', 'rs1000', 'ss2000'])

    @patch('tasks.eva_2357.populate_ids.get_mongo_uri_for_eva_profile')
    def test_populate_ids_fail(self, mock_get_mongo_uri_for_eva_profile):
        logging.getLogger().setLevel(logging.DEBUG)
        mock_get_mongo_uri_for_eva_profile.return_value = 'mongodb://127.0.0.1:27017'
        settings = self.get_test_resource("settings.xml")
        databases_file = self.get_test_resource("databases_errors.txt")
        with self.assertRaises(Exception):
            self.assertEqual(1, populate_ids(settings, databases_file, profile='localhost',
                                             mongo_accession_db=self.accession_db))

    @patch('tasks.eva_2357.populate_ids.get_mongo_uri_for_eva_profile')
    def test_check_all_contigs(self, mock_get_mongo_uri_for_eva_profile):
        logging.getLogger().setLevel(logging.DEBUG)
        mock_get_mongo_uri_for_eva_profile.return_value = 'mongodb://127.0.0.1:27017'
        settings = self.get_test_resource("settings.xml")
        databases_file = self.get_test_resource("databases.txt")
        check_all_contigs(settings, databases_file, profile='localhost')

    @patch('tasks.eva_2357.populate_ids.get_mongo_uri_for_eva_profile')
    def test_check_all_contigs_raise(self, mock_get_mongo_uri_for_eva_profile):
        logging.getLogger().setLevel(logging.DEBUG)
        mock_get_mongo_uri_for_eva_profile.return_value = 'mongodb://127.0.0.1:27017'
        settings = self.get_test_resource("settings.xml")
        databases_file = self.get_test_resource("databases_errors.txt")
        with pytest.raises(ValueError):
            check_all_contigs(settings, databases_file, profile='localhost')


    @staticmethod
    def get_test_resource(resource_name):
        """Gets full path to the test resource located in the same directory as the test module."""
        # Full path to this module.
        this_module = os.path.abspath(__file__)
        # Full path to the directory where it is contained.
        module_directory = os.path.dirname(this_module)
        # Full path to the requested resource.
        return os.path.join(module_directory, 'resources', resource_name)
