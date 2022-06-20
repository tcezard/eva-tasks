#!/usr/bin/env python
import os
import subprocess
from argparse import ArgumentParser

import yaml
from ebi_eva_common_pyutils import command_utils
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.config_utils import get_primary_mongo_creds_for_profile, get_accession_pg_creds_for_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query


logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def write_remapping_process_props_template(template_file_path):
    mongo_host, mongo_user, mongo_pass = get_primary_mongo_creds_for_profile(cfg['maven']['environment'],
                                                                             cfg['maven']['settings_file'])
    pg_url, pg_user, pg_pass = get_accession_pg_creds_for_profile(cfg['maven']['environment'],
                                                                  cfg['maven']['settings_file'])
    with open(template_file_path, 'w') as open_file:
        open_file.write(f'''spring.datasource.driver-class-name=org.postgresql.Driver
spring.datasource.url={pg_url}
spring.datasource.username={pg_user}
spring.datasource.password={pg_pass}
spring.datasource.tomcat.max-active=3

spring.jpa.generate-ddl=true

spring.data.mongodb.host={mongo_host}
spring.data.mongodb.port=27017
spring.data.mongodb.database=eva_accession_sharded
spring.data.mongodb.username={mongo_user}
spring.data.mongodb.password={mongo_pass}

spring.data.mongodb.authentication-database=admin
mongodb.read-preference=secondaryPreferred
spring.main.web-environment=false
spring.main.allow-bean-definition-overriding=true
spring.jpa.properties.hibernate.jdbc.lob.non_contextual_creation=true
logging.level.uk.ac.ebi.eva.accession.remapping=INFO
parameters.chunkSize=1000
''')
    return template_file_path


def process_one_taxonomy(assembly, taxid, scientific_name, target_assembly, resume):
    base_directory = cfg['remapping']['base_directory']
    nextflow_remapping_process = os.path.join(os.path.dirname(__file__), 'remediate.nf')
    assembly_directory = os.path.join(base_directory, taxid, assembly)
    work_dir = os.path.join(assembly_directory, 'work')
    prop_template_file = os.path.join(assembly_directory, 'template.properties')
    os.makedirs(work_dir, exist_ok=True)

    remapping_log = os.path.join(assembly_directory, 'remapping_process.log')
    remapping_config_file = os.path.join(assembly_directory, 'remapping_process_config_file.yaml')
    remapping_config = {
        'taxonomy_id': taxid,
        'source_assembly_accession': assembly,
        'target_assembly_accession': target_assembly,
        'species_name': scientific_name,
        'output_dir': assembly_directory,
        'genome_assembly_dir': cfg['genome_downloader']['output_directory'],
        'template_properties': write_remapping_process_props_template(prop_template_file),
        'remapping_config': cfg.config_file
    }
    for part in ['executable', 'nextflow', 'groovy']:
        remapping_config[part] = cfg[part]
    with open(remapping_config_file, 'w') as open_file:
        yaml.safe_dump(remapping_config, open_file)

    try:
        command = [
            cfg['executable']['nextflow'],
            '-log', remapping_log,
            'run', nextflow_remapping_process,
            '-params-file', remapping_config_file,
            '-work-dir', work_dir
        ]
        if resume:
            command.append('-resume')
        curr_working_dir = os.getcwd()
        os.chdir(assembly_directory)
        command_utils.run_command_with_output('Nextflow remapping process', ' '.join(command))
    except subprocess.CalledProcessError as e:
        logger.error('Nextflow remapping pipeline failed')
        raise e
    finally:
        os.chdir(curr_working_dir)


def process_one_assembly(assembly, resume):
    # Check the original remapping tracking table for the appropriate taxonomies and target assemblies
    query = (
        'SELECT taxonomy, scientific_name, assembly_accession'
        'FROM eva_progress_tracker.remapping_tracker '
        f"WHERE origin_assembly_accession='{assembly}'"
        "AND source='DBSNP'"
    )
    with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
        for taxid, scientific_name, target_assembly in get_all_results_for_query(pg_conn, query):
            process_one_taxonomy(assembly, taxid, scientific_name, target_assembly, resume)


def main():
    argparse = ArgumentParser(description='Run remediation for one assembly')
    argparse.add_argument('--assembly', help='Source assembly to remediate')
    argparse.add_argument('--resume', action='store_true', default=False)

    args = argparse.parse_args()
    cfg.load_config_file(os.getenv('REMAPPINGCONFIG'))
    
    process_one_assembly(args.assembly, args.resume)


if __name__ == "__main__":
    main()
