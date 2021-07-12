from unittest import TestCase

from ebi_eva_common_pyutils.mongodb import MongoDatabase

from tasks.eva_2480.swap_contig_with_correct_contig_in_study import swap_with_correct_contig


class TestSwapContigInStudy(TestCase):
    def setUp(self) -> None:
        self.db = "eva_accession_sharded_test"
        self.collection = "submittedVariantEntity"
        self.uri = "mongodb://localhost:27017/"
        self.mongo_source = MongoDatabase(uri=self.uri, db_name=self.db)
        contig_list = [
            {
                "_id": "D0D2E897BFD59EAF6E2D0BA2C7883B5DC5B34F30",
                "accession": 7169189340,
                "alt": "C",
                "contig": "AEMK02000229.1",
                "createdDate": "2020-09-01T22:55:50.956Z",
                "ref": "G",
                "seq": "GCA_000003025.6",
                "start": 24577,
                "study": "PRJEB28579",
                "tax": 9823,
                "version": 1
            },
            {
                "_id": "C5E53108D33B135EB45E2BFD3E67744B48CB06A8",
                "accession": 7166483872,
                "alt": "G",
                "contig": "AEMK02000626.1",
                "createdDate": "2020-09-01T15:31:27.489Z",
                "ref": "T",
                "seq": "GCA_000003025.6",
                "start": 15784,
                "study": "PRJEB28579",
                "tax": 9823,
                "version": 1
            }
        ]
        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle[self.db][self.collection].insert_many(contig_list)

    def tearDown(self) -> None:
        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle.close()

    def test_replace_with_correct_contig(self):
        contig_swap_list = [{'contig_1': 'AEMK02000229.1', 'contig_2': 'AEMK02000626.1'}]
        swap_with_correct_contig(self.mongo_source, contig_swap_list)

        variant1 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'study': 'PRJEB28579', 'seq': 'GCA_000003025.6', 'accession': 7169189340}))
        variant2 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'study': "PRJEB28579", 'seq': 'GCA_000003025.6', 'accession': 7166483872}))

        self.assertEqual(variant1['contig'], 'AEMK02000626.1')
        self.assertEqual(variant1['_id'], '35BCCD4FE99F1EFFA94704F0DA1141848468D361')
        self.assertEqual(variant2['contig'], 'AEMK02000229.1')
        self.assertEqual(variant2['_id'], '7859F0F6375E969174E4815B41D920BC601588C8')

        variant1 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': 'D0D2E897BFD59EAF6E2D0BA2C7883B5DC5B34F30'}))
        self.assertIsNone(variant1)
        variant2 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': 'C5E53108D33B135EB45E2BFD3E67744B48CB06A8'}))
        self.assertIsNone(variant2)
