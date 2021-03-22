# Copyright 2021 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import collections
import os
import yaml

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.nextflow import LinearNextFlowPipeline

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


class MoveMongoDBs:
    def __init__(self, migration_config_file, dbs_to_migrate_list, batch_number, resume_flag):
        self.migration_config = yaml.load(open(migration_config_file), Loader=yaml.FullLoader)
        self.dbs_to_migrate = [x.strip() for x in open(dbs_to_migrate_list).readlines() if x.strip() != '']
        self.batch_number = batch_number
        self.resume_flag = resume_flag
        # Processes, in order, that make up the workflow and the arguments that they take
        self.workflow_process_arguments_map = collections.OrderedDict(
            [("dump_data_from_source", ["mongo-source-uri", "mongo-source-secrets-file", "db-name", "dump-dir"]),
             ("prepare_dest_db", ["mongo-source-uri", "mongo-source-secrets-file",
                                  "mongo-dest-uri", "mongo-dest-secrets-file", "db-name"]),
             ("restore_data_to_dest", ["mongo-dest-uri", "mongo-dest-secrets-file", "db-name", "dump-dir"]),
             ("create_indexes_in_dest", ["mongo-source-uri", "mongo-source-secrets-file", "mongo-dest-uri",
                                         "mongo-dest-secrets-file", "db-name"])
             ])
        self.workflow_run_dir = os.path.join(self.migration_config["migration-folder"], f"batch{batch_number}")
        self.workflow_file_path = os.path.join(self.workflow_run_dir, f"batch{batch_number}_migration_workflow.nf")
        self.log_file_name = os.path.join(self.workflow_run_dir, f"batch{batch_number}_workflow_execution.log")
        os.makedirs(self.workflow_run_dir, exist_ok=True)

        self.db_move_pipeline = LinearNextFlowPipeline(workflow_file_path=self.workflow_file_path
                                                  , nextflow_binary_path=self.migration_config["nextflow-binary-path"]
                                                  , nextflow_config_path=self.migration_config["nextflow-config-path"]
                                                  , working_dir=self.workflow_run_dir)

    def add_processes_to_pipeline(self):
        for db_name in self.dbs_to_migrate:
            for workflow_process_name, workflow_process_args in self.workflow_process_arguments_map.items():
                process_config = {"dump-dir": self.workflow_run_dir, "db-name": db_name}
                process_with_args = "{0} {1}".format(workflow_process_name,
                                                     " ".join(["--{0} {1}".format(
                                                         arg,
                                                         self.migration_config.get(arg,
                                                                                   process_config.get(arg, '')))
                                                         for arg in workflow_process_args]))
                command_to_run = f"export PYTHONPATH={self.migration_config['script-path']} &&  " \
                                 f"({self.migration_config['python3-path']} " \
                                 f"-m eva_2338.{process_with_args} " \
                                 f"1>> {self.log_file_name} 2>&1)"
                self.db_move_pipeline.add_process(process_name=f"{workflow_process_name}_{db_name}",
                                                  command_to_run=command_to_run)

    def move(self):
        self.add_processes_to_pipeline()
        self.db_move_pipeline.run_pipeline(resume=self.resume_flag)


def main():
    parser = argparse.ArgumentParser(description='Move a database from one MongoDB instance to another',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--migration-config-file",
                        help="Full path to the migration configuration file (ex: /path/to/migration/config.yml)",
                        required=True)
    parser.add_argument("--dbs-to-migrate-list",
                        help="Full path to file with list of databases to migrate (ex: /path/to/dbs/to/migrate.txt)",
                        required=True)
    parser.add_argument("--batch", help="Migration batch (ex: 1)", required=True)
    parser.add_argument("--resume", help="Flag to indicate if migration job is to be resumed", action='store_true',
                        required=False)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()
    MoveMongoDBs(args.migration_config_file, args.dbs_to_migrate_list, args.batch, args.resume).move()


if __name__ == "__main__":
    main()
