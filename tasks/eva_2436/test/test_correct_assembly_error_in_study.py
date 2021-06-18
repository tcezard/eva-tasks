from unittest import TestCase

from ebi_eva_common_pyutils.mongodb import MongoDatabase

from tasks.eva_2436.correct_assembly_error_in_study import replace_with_correct_assembly


class TestCorrectContigErrorInStudy(TestCase):
    def setUp(self) -> None:
        self.assembly = "GCA_000002315.5"
        self.db = "eva_accession_sharded"
        self.collection = "submittedVariantEntity"
        self.uri = "mongodb://localhost:27017/"
        self.mongo_source = MongoDatabase(uri=self.uri, db_name=self.db)
        wrong_assembly = {
            "_id": "CF2DD190877BE29372CAB647EEC609525B861F22",
            "accession": 7054574707,
            "alt": "TTTGTCTGTGTATGGCTCTGGTGACACATCATTGCCTGGTGACAGGACTCT",
            "contig": "CM000094.5",
            "createdDate": "2020-01-27T10:53:40.169Z",
            "ref": "",
            "seq": "PRJEB36115",
            "start": 3672970,
            "study": "PRJEB36115",
            "tax": 9031,
            "version": 1
        }

        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle[self.db][self.collection].insert_one(wrong_assembly)

    def tearDown(self) -> None:
        self.mongo_source.mongo_handle[self.db][self.collection].drop()
        self.mongo_source.mongo_handle.close()

    def test_replace_with_correct_contig(self):
        fixed = replace_with_correct_assembly(self.mongo_source)
        self.assertEqual(fixed, 1)
        variant = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'study': 'PRJEB36115'}))
        self.assertEqual(variant['seq'], self.assembly)
        self.assertEqual(variant['_id'], '058F1637FC4817E8ABFDFFF3E501A0E1F5C41468')
        variant1 = (self.mongo_source.mongo_handle[self.db][self.collection].find_one(
            {'_id': 'CF2DD190877BE29372CAB647EEC609525B861F22'}))
        self.assertIsNone(variant1)
