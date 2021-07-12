from unittest import TestCase

from ebi_eva_common_pyutils.mongodb import MongoDatabase

from tasks.eva_2155.correct_contig_error_in_study import correct
import logging


class TestCorrectContigErrorInStudy(TestCase):
    def setUp(self) -> None:
        self.contig = "AF010406.1"
        self.db = "eva_accession_sharded_test"
        self.collection = "submittedVariantEntity"
        self.uri = "mongodb://localhost:27017/"
        self.mongo_source = MongoDatabase(uri=self.uri, db_name=self.db)
        wrong_contig = [{
                "_id": "8C4E490E82B895ADE3B9405204771B9E6BDCB286",
                "seq": "GCA_000298735.1",
                "tax": 9940,
                "study": "PRJEB33693",
                "contig": "-",
                "start": 16410,
                "ref": "G",
                "alt": "A",
                "accession": 7121896076,
                "version": 1,
                "createdDate": "2020-05-05T10:38:43.367Z"
            },
            {
                "_id": "7350E3A0B0242791BD25901F15D22467DC7939BD",
                "seq": "GCA_000298735.1",
                "tax": 9940,
                "study": "PRJEB33693",
                "contig": "OARMT",
                "start": 16410,
                "ref": "G",
                "alt": "A",
                "accession": 7121824383,
                "version": 1,
                "createdDate": "2020-04-28T00:26:01.844Z"
            },
            {
                "_id": "C9618202A2AF568A94259A1A16AB0A67DCC1CC94",
                "seq": "GCA_000298735.1",
                "tax": 9940,
                "study": "PRJEB23437",
                "contig": "CM001582.1",
                "start": 5442343,
                "ref": "C",
                "alt": "T",
                "accession": 5264373293,
                "version": 1,
                "createdDate": "2019-07-07T11:14:13.110Z"
            }

        ]

        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle[self.db][self.collection].insert_many(wrong_contig)

    def tearDown(self) -> None:
        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle.close()

    def test_correct(self):
        fixed = correct(self.mongo_source)
        self.assertEqual(fixed, (1, 1))
        updated_variant = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'seq': 'GCA_000298735.1', 'accession': 7121896076}))
        self.assertEqual(updated_variant['contig'], self.contig)
        self.assertEqual(updated_variant['_id'], '2FE89DBEFF0FB8CB4A544070042101916C9BA0A4')
        
        variant_with_hyphen = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': '8C4E490E82B895ADE3B9405204771B9E6BDCB286'}))
        self.assertIsNone(variant_with_hyphen)
        
        variant_with_OARMT = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': '7350E3A0B0242791BD25901F15D22467DC7939BD'}))
        self.assertIsNone(variant_with_OARMT)
        
        variant_to_keep = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': 'C9618202A2AF568A94259A1A16AB0A67DCC1CC94'}))
        self.assertIsNotNone(variant_to_keep)

        self.assertEqual(self.mongo_source.mongo_handle[self.db][self.collection].count(), 2)
