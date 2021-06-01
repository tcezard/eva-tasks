import json
import os.path
from collections import defaultdict
from itertools import islice

import psycopg2
import psycopg2.extras
from ebi_eva_common_pyutils.config_utils import get_pg_uri_for_variant_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongodb import MongoDatabase
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query

from migration_util import write_query_to_file

logger = logging_config.get_logger(__name__)

study_key = "input.study.id"
analysis_key = "input.vcf.id"
mongodb_key = 'spring.data.mongodb.database'
variant_collection = "variants_2_0"
files_collection = "files_2_0"
annotation_collection = "annotations_2_0"
annotation_metadata_collection = "annotationMetadata_2_0"
files_query_file_name = "files_query.txt"
variants_query_file_name = "variants_query.txt"
annotation_query_file_name = "annotations_query.txt"
annotation_metadata_query_file_name = "annotations_metadata_query.txt"
chunk_size = 1000


def find_variants_studies_eligible_for_migration(private_config_xml_file, migration_start_time, migration_end_time):
    with psycopg2.connect(get_pg_uri_for_variant_profile("production", private_config_xml_file),
                          user="evajt") as metadata_connection_handle:
        query_string = f"select bjep.job_execution_id, bjep.key_name, bjep.string_val, bje.start_time \
                        from batch_job_execution bje join batch_job_execution_params bjep \
                        on bje.job_execution_id=bjep.job_execution_id \
                        where bjep.key_name in ('{study_key}', '{analysis_key}', '{mongodb_key}')  \
                        and bje.start_time between '{migration_start_time}' and '{migration_end_time}'\
                        order by bjep.job_execution_id desc , bjep.key_name"

    query_result = get_all_results_for_query(metadata_connection_handle, query_string)
    logger.info(f"\nStudies eligible for migration : {query_result}")

    job_parameter_combine = defaultdict(dict)
    for job_id, key_name, key_value, start_time in query_result:
        job_parameter_combine[job_id][key_name] = key_value

    db_study_dict = defaultdict(set)
    for key, val in job_parameter_combine.items():
        db_study_dict[val[mongodb_key]].add((val[study_key], val[analysis_key]))

    return db_study_dict


def mongo_export_files_variants_data(mongo_source_uri, mongo_source_secrets_file, db_study_dict, export_dir, query_dir):
    logger.info(f"Starting mongo export process for  mongo ({mongo_source_uri})")
    for db, study_vcf in db_study_dict.items():
        mongo_source = MongoDatabase(uri=mongo_source_uri, secrets_file=mongo_source_secrets_file, db_name=db)
        files_query = create_files_query(study_vcf)
        files_query_path = write_query_to_file(files_query, query_dir, files_query_file_name)
        files_mongo_export_args = {
            "collection": files_collection,
            "queryFile": files_query_path
        }
        logger.info(
            f"Exporting data for database ({db}): collection ({files_collection}) - files_mongo_export_args ({files_mongo_export_args})")
        files_export_file = os.path.join(export_dir, db, files_collection, files_collection)
        mongo_source.export_data(files_export_file, files_mongo_export_args)

        variants_query = create_variants_query(study_vcf)
        variants_query_path = write_query_to_file(variants_query, query_dir, variants_query_file_name)
        variants_mongo_export_args = {
            "collection": variant_collection,
            "queryFile": variants_query_path
        }
        logger.info(
            f"Exporting data for database ({db}): collection ({variant_collection}) - variants_mongo_export_args ({variants_mongo_export_args})")
        variant_export_file = os.path.join(export_dir, db, variant_collection, variant_collection)
        mongo_source.export_data(variant_export_file, variants_mongo_export_args)


def annotations_export(mongo_source_uri, mongo_source_secrets_file, export_dir, query_dir):
    db_list = os.listdir(export_dir)
    for db in db_list:
        mongo_source = MongoDatabase(uri=mongo_source_uri, secrets_file=mongo_source_secrets_file, db_name=db)
        variant_file_loc = os.path.join(export_dir, db, variant_collection, variant_collection)
        if os.path.isfile(variant_file_loc):
            with open(variant_file_loc, 'r') as variant_file:
                chunk_number = 0
                while True:
                    variant_batch = list(islice(variant_file, chunk_size))
                    if not variant_batch:
                        break
                    annotations = get_annotations_ids(variant_batch)
                    annotation_ids = annotations["annotations_id"]
                    annotation_metadata_ids = annotations["annotations_metadata_id"]
                    if annotation_ids:
                        export_annotations_data(mongo_source, db, annotation_collection, annotation_ids, export_dir,
                                                query_dir, annotation_query_file_name, chunk_number)
                    if annotation_metadata_ids:
                        export_annotations_data(mongo_source, db, annotation_metadata_collection,
                                                annotation_metadata_ids, export_dir, query_dir,
                                                annotation_metadata_query_file_name, chunk_number)
                    chunk_number = chunk_number + 1


def export_annotations_data(mongo_source, db, collection, ids, export_dir, query_dir, query_file_name, chunk_number):
    query = create_query_with_ids(ids)
    query_file_path = write_query_to_file(query, query_dir, query_file_name)
    mongo_annot_export_args = {
        "collection": collection,
        "queryFile": query_file_path,
        "readPreference": "secondaryPreferred"
    }
    logger.info(
        f"Exporting data for database ({db} and collection ({collection}) - mongo_annot_export_args({mongo_annot_export_args})")
    export_file = os.path.join(export_dir, db, collection, f'{collection}_{chunk_number}')
    mongo_source.export_data(export_file, mongo_annot_export_args)


def get_annotations_ids(variant_batch):
    annotations_list = {
        "annotations_id": set(),
        "annotations_metadata_id": set()
    }
    for variant_str in variant_batch:
        variant = json.loads(variant_str)
        if "annot" not in variant:
            continue
        else:
            annot_array = variant["annot"]
            for annot in annot_array:
                annotations_list["annotations_id"].add(
                    json.dumps(f'{variant["_id"]}_{annot["vepv"]}_{annot["cachev"]}')[1:-1])
                annotations_list["annotations_metadata_id"].add(f'{annot["vepv"]}_{annot["cachev"]}')

    return annotations_list


def create_files_query(study_vcf):
    # query_template : {$or:[{sid:"PRJEB1234",fid:"EZFZV89"},{sid:"PRJEB4567",fid:"EVFZ067"},...]}
    query_beg = "{$or:["
    study_string = ",".join("{sid:\"" + x[0] + "\",fid:\"" + x[1] + "\"}" for x in study_vcf)
    query_end = "]}"
    files_query = query_beg + study_string + query_end
    logger.info(f"query created for files collection migration : {files_query}")

    return files_query


def create_variants_query(param_val):
    # query_template : {files:{$elemMatch:{\"$or\":[{sid:"PRJEB1234",fid:"EZFZV89"},{sid:"PRJEB4567",fid:"EVFZ067"},...]}}}
    query_beg = "{files:{$elemMatch:{\"$or\":["
    study_string = ",".join("{sid:\"" + x[0] + "\",fid:\"" + x[1] + "\"}" for x in param_val)
    query_end = "]}}}"
    variants_query = query_beg + study_string + query_end
    logger.info(f"query created for variants collection migration : {variants_query}")

    return variants_query


def create_query_with_ids(ids):
    # query_template : {_id:{$in:["1_4568_C_G_78_78","1_124578_A_T_86_86",....]}}
    query_beg = "{_id:{$in:["
    id_string = ",".join(f'"{x}"' for x in ids)
    query_end = "]}}"
    query_with_id = query_beg + id_string + query_end
    logger.info(f"query created for annotation migration : {query_with_id}")

    return query_with_id


def files_variants_export(mongo_source_uri, mongo_source_secrets_file, private_config_xml_file, export_dir,
                          query_file_dir, start_time, end_time):
    db_study_dict = find_variants_studies_eligible_for_migration(private_config_xml_file, start_time, end_time)
    mongo_export_files_variants_data(mongo_source_uri, mongo_source_secrets_file, db_study_dict, export_dir,
                                     query_file_dir)
