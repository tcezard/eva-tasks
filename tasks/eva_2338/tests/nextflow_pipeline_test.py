import os
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from tasks.eva_2338.move_mongo_dbs import MoveMongoDBs
from unittest import TestCase


class TestMoveMongoDBs(TestCase):
    # Tests require Nextflow binary and two locally running mongos instances in different ports
    def test_move(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        config_file_content = f"""migration-folder: {dir_path}/../resources
python3-path: python3
nextflow-binary-path: nextflow
nextflow-config-path: {dir_path}/workflow.config
script-path: {dir_path}/../../
mongo-source-uri: mongodb://localhost:27017/admin
mongo-source-secrets-file: {dir_path}/empty_secret_file 
mongo-dest-uri: mongodb://localhost:27018/admin
mongo-dest-secrets-file: {dir_path}/empty_secret_file
"""
        open(f"{dir_path}/migration_config.yml", "w").write(config_file_content)
        mover = MoveMongoDBs(migration_config_file=f"{dir_path}/migration_config.yml",
                             dbs_to_migrate_list=f"{dir_path}/dbs_to_migrate.txt",
                             batch_number="1", resume_flag=False)

        # Load data to source
        for db_name in mover.dbs_to_migrate:
            source_db = MongoDatabase(mover.migration_config["mongo-source-uri"], db_name=db_name)
            source_db.drop()
            source_db.restore_data(dump_dir=f"{dir_path}/../resources/{db_name}")

        mover.move()
        # Check if source data made it to the destination
        for db_name in mover.dbs_to_migrate:
            source_db = MongoDatabase(mover.migration_config["mongo-source-uri"], db_name=db_name)
            dest_db = MongoDatabase(mover.migration_config["mongo-dest-uri"], db_name=db_name)
            for collection_name in source_db.get_collection_names():
                self.assertEqual(source_db.mongo_handle[db_name][collection_name].count_documents(filter={}),
                                 dest_db.mongo_handle[db_name][collection_name].count_documents(filter={}))
