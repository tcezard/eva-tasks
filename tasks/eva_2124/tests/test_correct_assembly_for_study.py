from datetime import datetime
from unittest import TestCase

import pymongo
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle

from tasks.eva_2124.correct_assembly_for_study import correct, get_mongo_connection_handle_url


class TestCorrectAssembly(TestCase):

    def setUp(self) -> None:
        # setup mongo database
        self.host = 'localhost'
        self.port = 27017
        self.eva_database = 'eva_accession_sharded'
        self.variant_collection = 'submittedVariantEntity'
        self.connection_handle = get_mongo_connection_handle_url(self.host, self.port)
        # connection_handle[self.admin_db].command("dropUser", self.username)
        # connection_handle[self.admin_db].command("createUser", self.username, pwd=self.password,
        #                                          roles=["userAdminAnyDatabase"])
        # connection_handle.close()

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
        horse_doc = {
            '_id': '0000185B0438E7CB2DD5E8F9F2B62E89C398FDA4',
            'seq': 'GCA_000002305.1',
            'tax': 9796,
            'study': 'PRJEB9799',
            'contig': '1',
            'start': 4573797,
            'ref': 'G',
            'alt': 'A',
            'accession': 5000016398,
            'version': 1,
            'createdDate': datetime.fromisoformat('2018-06-20T08:31:07.900')
        }
        another_horse_doc = {
            '_id' : 'D81824774C9D92D71971A970567FF9C5F510D283',
            'seq' : 'GCA_000002305.1',
            'tax' : 9796,
            'study' : 'PRJEB27771',
            'contig' : '21',
            'start' : 30681601,
            'ref' : 'C',
            'alt' : 'T',
            'accession' : 5081976263,
            'version' : 1,
            'createdDate' : datetime.fromisoformat('2018-09-20T14:31:26.164')
        }

        horse_already_genbank = {
            '_id' : 'BB6C2E844664691DC42FC46456448833A8A35763',
            'seq' : 'GCA_000002305.1',
            'tax' : 9796,
            'study' : 'PRJEB24630',
            'contig' : 'CM000387.2',
            'start' : 15520392,
            'ref' : 'C',
            'alt' : 'T',
            'accession' : 5201460111,
            'version' : 1,
            'createdDate' : datetime.fromisoformat('2019-07-08T03:12:59.352')
        }
        # self.connection_handle = get_mongo_connection_handle(username=self.username, password=self.password,
        #                                                      host=self.host, authentication_database=self.admin_db)
        self.connection_handle[self.eva_database][self.variant_collection].drop()
        self.connection_handle[self.eva_database][self.variant_collection].insert_many([
            original_document, horse_doc, another_horse_doc, horse_already_genbank])

    def tearDown(self) -> None:
        self.connection_handle[self.eva_database][self.variant_collection].drop()
        # self.connection_handle[self.admin_db].command("dropUser", self.username)
        self.connection_handle.close()

    def test_correct(self):
        fixed = correct(None, None, self.host, studies=['PRJEB30080'],
                        assembly_accession='GCA_000181335.3')
        self.assertEqual(fixed, 1)
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one({'seq': 'GCA_000181335.3'}))
        self.assertEqual(variant['contig'], 'CM001381.2')
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one(
            {'seq': 'GCA_000181335.3', 'contig': 'B1'}))
        self.assertIsNone(variant)

    def test_already_genbank(self):
        fixed = correct(None, None, self.host, studies=['PRJEB9799', 'PRJEB24630', 'PRJEB27771'],
                        assembly_accession='GCA_000002305.1')
        self.assertEqual(fixed, 2)
        variants = (self.connection_handle[self.eva_database][self.variant_collection].find({'seq': 'GCA_000002305.1'}))
        for variant in variants:
            self.assertTrue(variant['contig'] in ['CM000377.2', 'CM000397.2', 'CM000387.2'])
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one(
            {'seq': 'GCA_000181335.3', 'contig': '1'}))
        self.assertIsNone(variant)

    def test_ignore_otherstudies(self):
        fixed = correct(None, None, self.host, studies=['PRJEB9799'],
                        assembly_accession='GCA_000002305.1')
        self.assertEqual(fixed, 1)
        variants = (self.connection_handle[self.eva_database][self.variant_collection].find({'seq': 'GCA_000002305.1'}))
        for variant in variants:
            self.assertTrue(variant['contig'] in ['CM000377.2', '21', 'CM000387.2'])

    def test_chunk_size(self):
        fixed = correct(None, None, self.host, studies=['PRJEB9799', 'PRJEB24630', 'PRJEB27771'],
                        assembly_accession='GCA_000002305.1', chunk_size=1)
        self.assertEqual(fixed, 2)
        variants = (self.connection_handle[self.eva_database][self.variant_collection].find({'seq': 'GCA_000002305.1'}))
        for variant in variants:
            self.assertTrue(variant['contig'] in ['CM000377.2', 'CM000397.2', 'CM000387.2'])
        variant = (self.connection_handle[self.eva_database][self.variant_collection].find_one(
            {'seq': 'GCA_000181335.3', 'contig': '1'}))
        self.assertIsNone(variant)

    def test_can_not_replace_contig(self):
        unreplaceable_doc = {
            '_id': 'B0B2176A78EDD8B29FEC60E87941287DC9E5056B',
            'seq': 'GCA_000002305.1',
            'tax': 9796,
            'study': 'PRJEB27771',
            'contig': 'asdf',
            'start': 30681601,
            'ref': 'C',
            'alt': 'T',
            'accession': 5081976263,
            'version': 1,
            'createdDate': datetime.fromisoformat('2018-09-20T14:31:26.164')
        }
        self.connection_handle[self.eva_database][self.variant_collection].insert_one(unreplaceable_doc)
        try:
            fixed = correct(None, None, self.host, studies=['PRJEB9799', 'PRJEB24630', 'PRJEB27771'],
                            assembly_accession='GCA_000002305.1')
            self.fail("expected an exception but none was thrown")
        except Exception:
            variants = (self.connection_handle[self.eva_database][self.variant_collection].find(
                {'seq': 'GCA_000002305.1'}))
            original_contigs_before_replacement = ['CM000387.2', '1', '21', 'asdf']
            for variant in variants:
                self.assertIn(variant['contig'], original_contigs_before_replacement)
