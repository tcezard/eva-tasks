from mock import patch

from pymongo import MongoClient
from unittest import TestCase
from tasks.eva_2287.correct_contig_with_chr_removed import correct


class TestCorrectChr(TestCase):
    def setUp(self) -> None:
        self.host = 'localhost'
        self.port = 27017
        self.database = 'GCA_002742125_1_test'
        self.variant_collection = 'submittedVariantEntity'
        self.connection_handle = MongoClient(self.host)

        wrong_contig = {
            "_id": "616766226DF5A3852B9FE4B266B641421CAF4BE6",
            "seq": "GCA_002742125.1",
            "tax": 9940,
            "study": "PRJEB42582",
            "contig": "omosome11",
            "start": 25929637,
            "ref": "G",
            "alt": "A",
            "accession": 7180281342,
            "version": 1,
            "createdDate": "2021-01-26T16:56:08.555Z"
        }

        self.connection_handle[self.database][self.variant_collection].drop()
        self.connection_handle[self.database][self.variant_collection].insert_many([wrong_contig])

    def tearDown(self) -> None:
        self.connection_handle[self.database][self.variant_collection].drop()
        self.connection_handle.close()

    @patch('tasks.eva_2287.correct_contig_with_chr_removed.get_mongo_uri_for_eva_profile')
    def test_correct(self, mock_get_mongo_uri_for_eva_profile):
        mock_get_mongo_uri_for_eva_profile.return_value = 'mongodb://127.0.0.1:27017'
        fixed = correct('../test/settings.xml', profile='localhost', mongo_database=self.database)
        self.assertEqual(fixed, 1)
        variant = (self.connection_handle[self.database][self.variant_collection].find_one(
            {'seq': 'GCA_002742125.1'}))
        self.assertEqual(variant['contig'], 'chromosome11')
        self.assertEqual(variant['_id'], '55392EB026F670F15823A9BEF16DD45398371124')
        variant = (self.connection_handle[self.database][self.variant_collection].find_one(
            {'_id': '616766226DF5A3852B9FE4B266B641421CAF4BE6'}))
        self.assertIsNone(variant)
