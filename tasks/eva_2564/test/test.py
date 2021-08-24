from unittest import TestCase

from ebi_eva_common_pyutils.mongodb import MongoDatabase

from tasks.eva_2564.correct_contig_error_in_study import correct


class TestCorrectContigErrorInStudy(TestCase):
    def setUp(self) -> None:
        self.db = "eva_accession_sharded_test"
        self.collection = "submittedVariantEntity"
        self.uri = "mongodb://localhost:27017/"
        self.mongo_source = MongoDatabase(uri=self.uri, db_name=self.db)
        wrong_contig = [
            {
                "_id": "52AF357592EAE966B65D3F9D04C952791F1493DC",
                "seq": "GCA_000001895.4",
                "tax": 10116,
                "study": "PRJEB42012",
                "contig": "Un.1",
                "start": 1063,
                "ref": "G",
                "alt": "C",
                "accession": 7315398622,
                "version": 1,
                "createdDate": "2021-02-09T03:23:05.842Z"
            },
            {
                "_id": "899D3B4E7B6E1B18B9D34AC0CA3880AF5F68E09B",
                "seq": "GCA_000001895.4",
                "tax": 10116,
                "study": "PRJEB42012",
                "contig": "1_random.1",
                "start": 1510,
                "ref": "T",
                "alt": "C",
                "accession": 7315398494,
                "version": 1,
                "createdDate": "2021-02-09T03:23:05.041Z"
            }
        ]

        correct_contig = [
            {
                "_id": "CAB0D97D36233AC8D84637F228B1CE172228A166",
                "seq": "GCA_000001895.4",
                "tax": 10116,
                "study": "PRJEB42012",
                "contig": "CM000072.5",
                "start": 2203542,
                "ref": "T",
                "alt": "C",
                "accession": 7306435312,
                "version": 1,
                "createdDate": "2021-02-08T11:56:42.206Z"
            }
        ]

        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle[self.db][self.collection].insert_many(wrong_contig)
        self.mongo_source.mongo_handle[self.db][self.collection].insert_many(correct_contig)

    def tearDown(self) -> None:
        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle.close()

    def test_correct(self):
        fixed = correct(self.mongo_source)
        self.assertEqual(fixed, (2, 2))

        updated_variant_1 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'seq': 'GCA_000001895.4', 'accession': 7315398622}))
        self.assertEqual(updated_variant_1['contig'], 'KL568167.1')
        self.assertEqual(updated_variant_1['_id'], '5633D310C37E5FE1AA374A9D203CE91471ADCF02')

        updated_variant_2 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'seq': 'GCA_000001895.4', 'accession': 7315398494}))
        self.assertEqual(updated_variant_2['contig'], 'AABR07046142.1')
        self.assertEqual(updated_variant_2['_id'], '48E11C28F20EFBD8DE8AF730AF523A2C71322730')

        variant_wrong_1 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': '52AF357592EAE966B65D3F9D04C952791F1493DC'}))
        self.assertIsNone(variant_wrong_1)

        variant_wrong_2 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': '899D3B4E7B6E1B18B9D34AC0CA3880AF5F68E09B'}))
        self.assertIsNone(variant_wrong_2)

        variant_to_keep = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': '5633D310C37E5FE1AA374A9D203CE91471ADCF02'}))
        self.assertIsNotNone(variant_to_keep)

        self.assertEqual(self.mongo_source.mongo_handle[self.db][self.collection].count_documents({}), 3)
