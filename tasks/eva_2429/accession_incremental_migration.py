import os.path
from collections import defaultdict

import psycopg2
import psycopg2.extras
from ebi_eva_common_pyutils.config_utils import get_pg_uri_for_accession_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query

from migration_util import write_query_to_file

logger = logging_config.get_logger(__name__)

project_key = "projectAccession"
assembly_key = "assemblyAccession"
accession_query_file_name = "accession_query.txt"
accession_db = "eva_accession_sharded"
accession_collection = "submittedVariantEntity"


def find_accession_studies_eligible_for_migration(private_config_xml_file, migration_start_time, migration_end_time):
    with psycopg2.connect(get_pg_uri_for_accession_profile("production", private_config_xml_file),
                          user="evaaccjt") as metadata_connection_handle:
        query_string = f"select bjep.job_execution_id, bjep.key_name, bjep.string_val, bje.start_time \
                     from batch_job_execution bje join batch_job_execution_params bjep \
                     on bje.job_execution_id=bjep.job_execution_id \
                     where bjep.key_name in ('{project_key}', '{assembly_key}')  \
                     and bje.start_time between '{migration_start_time}' and '{migration_end_time}' \
                     order by bjep.job_execution_id desc , bjep.key_name"

    query_result = get_all_results_for_query(metadata_connection_handle, query_string)
    logger.info(f"\nStudies eligible for migration : {query_result}")

    job_parameter_combine = defaultdict(dict)
    for job_id, key_name, key_value, start_time in query_result:
        job_parameter_combine[job_id][key_name] = key_value

    study_seq_tuple_set = set()
    for key, val in job_parameter_combine.items():
        study_seq_tuple_set.add((val[project_key], val[assembly_key]))

    return study_seq_tuple_set


def export_accession_data(mongo_source_uri, mongo_source_secrets_file, study_seq_tuple_set, export_dir, query_file_dir):
    mongo_source = MongoDatabase(uri=mongo_source_uri, secrets_file=mongo_source_secrets_file, db_name=accession_db)
    accession_query = create_accession_query(study_seq_tuple_set)
    query_file_path = write_query_to_file(accession_query, query_file_dir, accession_query_file_name)
    mongo_export_args = {
        "collection": accession_collection,
        "queryFile": query_file_path
    }

    logger.info(
        f"Starting mongo export process for accessioning database: mongo_source ({mongo_source_uri}) and mongo_export_args ({mongo_export_args})")
    accession_export_file = os.path.join(export_dir, accession_db, accession_collection, accession_collection)

    mongo_source.export_data(accession_export_file, mongo_export_args)


def create_accession_query(study_seq_tuple_set):
    # query_template: {$or:[{study:"PRJEB1234",seq:"GC000012.4"},{study:"PRJEB4567",seq:"GC000045.3"},...]}
    query_beg = "{$or:["
    study_string = ",".join("{study:\"" + x[0] + "\",seq:\"" + x[1] + "\"}" for x in study_seq_tuple_set)
    query_end = "]}"
    accession_query = query_beg + study_string + query_end
    logger.info(f"query created for accession migration : {accession_query}")

    return accession_query


def accession_export(mongo_source_uri, mongo_source_secrets_file, private_config_xml_file, export_dir, query_file_dir,
                     start_time, end_time):
    study_seq_set = find_accession_studies_eligible_for_migration(private_config_xml_file, start_time, end_time)
    export_accession_data(mongo_source_uri, mongo_source_secrets_file, study_seq_set, export_dir, query_file_dir)
