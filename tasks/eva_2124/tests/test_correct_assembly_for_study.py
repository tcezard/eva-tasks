from datetime import datetime
from unittest import TestCase

import pymongo
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle

from tasks.eva_2124.correct_assembly_for_study import correct


class TestCorrectAssembly(TestCase):

    def setUp(self) -> None:
        # setup mongo database
        self.host = 'localhost'
        self.port = 27017
        self.admin_db = 'admin'
        self.username = 'admin'
        self.password = 'admin'
        self.eva_database = 'eva_accession_sharded'
        self.variant_collection = 'submittedVariantEntity'
        mongo_connection_uri = "mongodb://{0}:{1}/{2}".format(self.host, self.port, self.admin_db)
        connection_handle = pymongo.MongoClient(mongo_connection_uri)
        connection_handle[self.admin_db].command("createUser", self.username, pwd=self.password,
                                                 roles=["userAdminAnyDatabase"])
        connection_handle.close()

        original_document = {
            '_id': '16261B7E16FD9AC7CB79D640736000F3C9FE2122',
            'seq': 'GCA_000181335.3',
            'tax': 9685,
            'study': 'PRJEB30080',
            'contig': 'B1',
            'start': 150039652,
            'ref': 'C',
            'alt': 'T',
            'accession': 5015497292,
            'version': 1,
            'createdDate': datetime.fromisoformat('2018-12-05T14:18:37.842')
        }

        self.connection_handle = get_mongo_connection_handle(username=self.username, password=self.password,
                                                             host=self.host, authentication_database=self.admin_db)
        self.connection_handle[self.eva_database][self.variant_collection].insert_one(original_document)

    def tearDown(self) -> None:
        self.connection_handle[self.eva_database][self.variant_collection].drop()
        self.connection_handle[self.admin_db].command("dropUser", self.username)
        self.connection_handle.close()

    def test_correct(self):
        correct(self.username, self.password, self.host, study='PRJEB33514', reference_source='GCA_000002305.1')
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one({'seq': 'GCA_000181335.3'}))
        self.assertEqual(variant['contig'], 'CM001381.2')
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one(
            {'seq': 'GCA_000181335.3', 'contig': 'B1'}))
        self.assertIsNone(variant)

