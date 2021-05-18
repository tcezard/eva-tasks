import argparse
from datetime import datetime

from ebi_eva_common_pyutils.logger import logging_config

from accession_incremental_migration import accession_export
from migration_util import mongo_import_from_dir
from variant_incremental_migration import variants_export

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()

all_tasks = ['accession_export', 'variant_export', 'import']


def main():
    parser = argparse.ArgumentParser(
        description='Given a mongo instance and a timestamp find the documents inserted in '
                    'variant databases after the timestamp and migrate them',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--mongo-source-secrets-file",
                        help="Full path to the Mongo Source secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--export-dir", help="Top-level directory where all export reside (ex: /path/to/archives)",
                        required=True)
    parser.add_argument("--mongo-dest-uri",
                        help="Mongo Destination URI (ex: mongodb://user:@mongos-source-host:27017/admin)",
                        required=True)
    parser.add_argument("--mongo-dest-secrets-file",
                        help="Full path to the Mongo Destination secrets file (ex: /path/to/mongo/source/secret)",
                        required=True)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml",
                        required=True)
    parser.add_argument("--start-time",
                        help="Studies that were processed after start time will be migrated (ex: \"2021-04-06 12:00:00\")",
                        required=True)
    parser.add_argument("--end-time",
                        help="Studies that were processed before end time will be migrated(ex: \"2021-04-06 12:00:00.000000\") \
                             If not provided current timestamp will be taken by default",
                        required=False, default=str(datetime.now()))
    parser.add_argument('--tasks', required=False, type=str, nargs='+',
                        default=all_tasks, choices=all_tasks,
                        help='Task or set of tasks to perform during migration.')
    parser.add_argument("--query-file-dir",
                        help="Top level directory where all the query files will be created. If not provided, script directory path will be taked by default",
                        required=False)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    if 'accession_export' in args.tasks:
        accession_export(args.mongo_source_uri, args.mongo_source_secrets_file, args.private_config_xml_file,
                         args.export_dir, args.query_file_dir, args.start_time, args.end_time)
    if 'variant_export' in args.tasks:
        variants_export(args.mongo_source_uri, args.mongo_source_secrets_file, args.private_config_xml_file,
                        args.export_dir, args.query_file_dir, args.start_time, args.end_time)
    if 'import' in args.tasks:
        mongo_import_from_dir(args.mongo_dest_uri, args.mongo_dest_secrets_file, args.export_dir)


if __name__ == "__main__":
    main()
