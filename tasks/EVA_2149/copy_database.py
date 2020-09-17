from argparse import ArgumentParser

from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
from ebi_eva_common_pyutils.mongo_utils import copy_db
from pymongo.uri_parser import parse_uri

logger = log_cfg.get_logger(__name__)


def copy_database_to(source_database, destination_database, private_config_xml_file, dump_dir):
    mongo_params = parse_uri(get_mongo_uri_for_eva_profile("production", private_config_xml_file))
    # nodelist is in format: [(host1,port1), (host2,port2)]. Just choose one.
    # Mongo is smart enough to fallback to secondaries automatically.
    mongo_host = mongo_params["nodelist"][0][0]
    logger.info("Beginning data copy for: " + source_database)
    dump_output_dir = "{0}/dump_{1}".format(dump_dir, source_database.replace(".", "_"))

    mongodump_args = {"db": source_database, "host": mongo_host,
                      "username": mongo_params["username"], "password": mongo_params["password"],
                      "authenticationDatabase": "admin", "out": dump_output_dir
                      }
    mongorestore_args = {"db": destination_database, "host": mongo_host,
                         "username": mongo_params["username"], "password": mongo_params["password"],
                         "dir": "{0}/{1}/".format(dump_output_dir, source_database)}
    logger.info("Running export to {0} against {2} in production".format(dump_output_dir, source_database ))
    copy_db(mongodump_args, mongorestore_args)


if __name__ == "__main__":
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--source_database')
    arg_parser.add_argument('--destination_database')
    arg_parser.add_argument('--dump_dir')
    args = arg_parser.parse_args()

    copy_database_to(args.source_database, args.destination_database, args.dump_dir)
