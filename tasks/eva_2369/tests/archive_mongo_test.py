import os
import sys
import tempfile
from unittest import TestCase

import pymongo
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.mongodb import MongoDatabase

from tasks.eva_2369.archive_mongo_dbs import main


class TestArchiveMongoDatabase(TestCase):
    db_names_list_file = "db-names-list.txt"
    mongo_secret_file = "mongo-secret-file.txt"
    db_name = "test_mongo_db"
    uri = "mongodb://localhost:27017/admin"
    dir_path = os.path.dirname(os.path.realpath(__file__))
    resources_folder = f"""{dir_path}/../resources"""

    local_mongo_handle = pymongo.MongoClient()

    # Tests expect a local sharded Mongo instance
    def setUp(self) -> None:
        self.test_mongo_db = MongoDatabase(uri=self.uri, db_name=self.db_name)
        self.dump_dir = os.path.join(self.resources_folder, self.db_name)
        run_command_with_output("Drop target test database if it already exists...", f"mongo {self.db_name} "
                                                                                     f"--eval 'db.dropDatabase()'")
        run_command_with_output("Import test database...", f"mongorestore --dir {self.dump_dir}")

    def test_archive_data(self):
        with tempfile.TemporaryDirectory() as tempdir:
            sys.argv[1:] = [f"""--db-names-list-file={self.dir_path}/../resources/{self.db_names_list_file}""",
                            f"""--archive-dir={tempdir}""",
                            f"""--mongo-source-uri={self.uri}""",
                            f"""--mongo-source-secrets-file={self.dir_path}/../resources/{self.mongo_secret_file}"""]
            main()
            self.assertTrue(os.path.isfile(os.path.join(tempdir, self.db_name)))
