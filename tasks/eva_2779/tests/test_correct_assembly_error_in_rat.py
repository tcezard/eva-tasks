import hashlib

from ebi_eva_common_pyutils.mongodb import MongoDatabase
from unittest import TestCase

from tasks.eva_2779.correct_assembly_error_in_rat import get_submitted_SHA1, get_clustered_SHA1, \
    replace_document_with_correct_information, update_operation_entities


def calculate_id(rs):
    h = hashlib.sha1()
    h.update(str(rs).encode())
    return h.hexdigest().upper()


class TestCorrectAssemblyError(TestCase):
    def setUp(self) -> None:
        host = 'localhost'
        port = 27017

        # Accessioning warehouse
        self.accession_db = 'eva_accession_sharded'
        self.submitted_variants_collection = 'submittedVariantEntity'
        self.clustered_variants_collection = 'clusteredVariantEntity'
        self.clustered_variants_operation_collection = 'clusteredVariantOperationEntity'
        self.dbsnp_clustered_variants_collection = 'dbsnpClusteredVariantEntity'
        uri = f'mongodb://{host}:{port}'
        self.mongo_db = MongoDatabase(uri=uri, db_name="eva_accession_sharded")
        self.connection_handle = self.mongo_db.mongo_handle

        # Accessioning db
        ss1 = {
            "seq": "GCA_015227675.1",
            "tax": 1111,
            "study": "PRJEB30318",
            "contig": "CM026996.1",
            "start": 76166296,
            "ref": "C",
            "alt": "T",
            "accession": 5318166021,
            "rs": 1000,
        }
        ss1["_id"] = get_submitted_SHA1(ss1)

        rs1 = {
            "asm": "GCA_015227675.1",
            "tax": 1111,
            "contig": "CM026996.1",
            "start": 76166296,
            "accession": 1000,
            'type': 'SNV'
        }
        rs1["_id"] = get_clustered_SHA1(rs1)

        cvoe1 = {
            "_id": "61a8fe72b0cc848d7a5a194f",
            "eventType": "MERGED",
            "accession": 198834536,
            "mergeInto": 198362526,
            "reason": "After remapping to GCA_015227675.1, RS IDs mapped to the same locus.",
            "inactiveObjects": [
                {
                        "asm": "GCA_015227675.1",
                        "tax": 10116,
                        "contig": "CM026996.1",
                        "start": 39058078,
                        "type": "SNV",
                        "hashedMessage": "AF611B8ED3BE813B427AD7E48D266D3DB5F013B8",
                        "accession": 198834536,
                        "version": 1,
                },
                {
                    "asm": "GCA_015227675.1",
                    "tax": 10116,
                    "contig": "CM026996.1",
                    "start": 39058077,
                    "type": "SNV",
                    "hashedMessage": "AF611B8ED3BE813B427AD7E48D266D3DB5F013B8",
                    "accession": 198834532,
                    "version": 1,
                }
            ]
        }

        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle[self.accession_db][self.clustered_variants_collection].drop()
        self.connection_handle[self.accession_db][self.clustered_variants_operation_collection].drop()
        self.connection_handle[self.accession_db][self.submitted_variants_collection].insert_many([ss1])
        self.connection_handle[self.accession_db][self.clustered_variants_collection].insert_many([rs1])
        self.connection_handle[self.accession_db][self.clustered_variants_operation_collection].insert_many([cvoe1])

    def tearDown(self) -> None:
        self.connection_handle[self.accession_db][self.submitted_variants_collection].drop()
        self.connection_handle[self.accession_db][self.clustered_variants_collection].drop()
        self.connection_handle.close()

    def test_replace_document_with_correct_information(self):
        replace_document_with_correct_information(
            self.mongo_db, "submittedVariantEntity", get_submitted_SHA1,
            {'seq': 'GCA_015227675.1', 'contig': 'CM026996.1'},
            {'seq': 'GCA_015227675.2', 'contig': 'AY172581.1'},
            5
        )
        cursor = self.connection_handle[self.accession_db][self.submitted_variants_collection].find()
        document = cursor.__next__()
        assert document['seq'] == 'GCA_015227675.2'
        assert document['contig'] == 'AY172581.1'

    def test_update_operation_entities(self):
        update_operation_entities(
            self.mongo_db, 'clusteredVariantOperationEntity', get_clustered_SHA1,
            {'inactiveObjects.asm': 'GCA_015227675.1', 'eventType': 'MERGED'},
            {'inactiveObjects.asm': 'GCA_015227675.2',
             'reason': lambda x: x.replace('GCA_015227675.1', 'GCA_015227675.2')},
            5
        )
        cursor = self.connection_handle[self.accession_db][self.clustered_variants_operation_collection].find()
        document = cursor.__next__()
        assert document['reason'] == 'After remapping to GCA_015227675.2, RS IDs mapped to the same locus.'
        assert document['inactiveObjects'][0]['asm'] == 'GCA_015227675.2'
