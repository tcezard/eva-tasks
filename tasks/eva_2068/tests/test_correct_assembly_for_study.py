from datetime import datetime
from unittest import TestCase

import pymongo
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle

from tasks.eva_2068.correct_assembly_for_study import correct


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
            '_id': '098B7F0FD7D5D99FD9B1514FBE4B5C5A8294D561',
            'seq': 'GCA_000002325.1',
            'tax': 7425,
            'study': 'PRJEB33514',
            'contig': 'CM000915.2',
            'start': 51524,
            'ref': 'T',
            'alt': 'A',
            'accession': 7123539069,
            'version': 1,
            'createdDate': datetime.fromisoformat('2020-06-24T15:20:46.533')
        }

        self.connection_handle = get_mongo_connection_handle(username=self.username, password=self.password,
                                                             host=self.host, authentication_database=self.admin_db)
        self.connection_handle[self.eva_database][self.variant_collection].insert_one(original_document)

    def tearDown(self) -> None:
        self.connection_handle[self.eva_database][self.variant_collection].drop()
        self.connection_handle[self.admin_db].command("dropUser", self.username)
        self.connection_handle.close()

    def test_correct(self):
        correct(self.username, self.password, self.host,
                study='PRJEB33514', reference_source='GCA_000002325.1',
                reference_dest='GCA_000002325.2')
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one({'seq': 'GCA_000002325.2'}))
        self.assertEqual(variant['_id'], 'CB5637C29FD1A9C3FA64D0C1EE2FC8B9147DEC47')
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one({'seq': 'GCA_000002325.1'}))
        self.assertIsNone(variant)

