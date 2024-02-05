import argparse
from urllib.parse import urlsplit

import psycopg2
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
from ebi_eva_internal_pyutils.config_utils import get_properties_from_xml_file
from ebi_eva_internal_pyutils.pg_utils import get_all_results_for_query, execute_query, logger

logger = log_cfg.get_logger(__name__)


def get_contig_alias_db_cred_for_profile(settings_xml_file, profile):
    properties = get_properties_from_xml_file(profile, settings_xml_file)
    contig_alias_db = properties['contig-alias.db-url']
    contig_alias_user = properties['contig-alias.db-username']
    contig_alias_pass = properties['contig-alias.db-password']
    return contig_alias_db, contig_alias_user, contig_alias_pass


def get_contig_alias_connection_handle(settings_xml_file, profile):
    pg_url, pg_user, pg_pass = get_contig_alias_db_cred_for_profile(settings_xml_file, profile)
    return psycopg2.connect(urlsplit(pg_url).path, user=pg_user, password=pg_pass)


def get_assemblies_to_update(private_config_xml_file, source_env, target_env, assembly_list):
    target_asm = set()
    target_query = f"select distinct assembly_insdc_accession from chromosome"
    with get_contig_alias_connection_handle(private_config_xml_file, target_env) as target_db_conn:
        for assembly in get_all_results_for_query(target_db_conn, target_query):
            target_asm.add(assembly[0])

    if assembly_list:
        assemblies = assembly_list[0].split(',')
        common_asm = target_asm.intersection(assemblies)
        return [asm for asm in common_asm]
    else:
        src_asm = []
        src_query = f"""select assembly_insdc_accession, count(*) as count from chromosome 
                    where md5checksum is not null 
                    and assembly_insdc_accession in {"(" + ",".join(["'" + asm + "'" for asm in target_asm]) + ")"} 
                    group by assembly_insdc_accession order by count"""
        with get_contig_alias_connection_handle(private_config_xml_file, source_env) as source_db_conn:
            for assembly, count in get_all_results_for_query(source_db_conn, src_query):
                src_asm.append(assembly)

        return src_asm


def copy_md5checksum_for_assemblies(private_config_xml_file, source_env, target_env, assemblies):
    with get_contig_alias_connection_handle(private_config_xml_file, source_env) as source_db_conn:
        with get_contig_alias_connection_handle(private_config_xml_file, target_env) as target_db_conn:
            for asm in assemblies:
                logger.info(f"Start updating MD5 checksum for assembly: {asm}")
                src_query = f"""select insdc_accession, md5checksum from chromosome 
                                where md5checksum is not null and assembly_insdc_accession='{asm}'"""
                for insdc_acc, md5checksum in get_all_results_for_query(source_db_conn, src_query):
                    target_query = f"""update chromosome set md5checksum='{md5checksum}' 
                                        where assembly_insdc_accession='{asm}' and insdc_accession='{insdc_acc}'"""
                    execute_query(target_db_conn, target_query)


def copy_md5checksum_from_source_to_prod(private_config_xml_file, source_env, target_env, assembly_list):
    assemblies = get_assemblies_to_update(private_config_xml_file, source_env, target_env, assembly_list)
    logger.info(f"Updating MD5 checksum for assemblies: {assemblies}")
    copy_md5checksum_for_assemblies(private_config_xml_file, source_env, target_env, assemblies)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Copy MD5 checksum from one env to another', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--source-env", choices=('localhost', 'development', 'production'),
                        help="Source env to copy data from", required=True)
    parser.add_argument("--target-env", choices=('localhost', 'development', 'production'),
                        help="Target env to copy data to", required=True)
    parser.add_argument("--assembly-list", help="Comma separated assembly list e.g. GCA_000181335.4,GCA_000181335.5",
                        required=False, nargs='+')

    args = parser.parse_args()

    copy_md5checksum_from_source_to_prod(args.private_config_xml_file, args.source_env, args.target_env,
                                         args.assembly_list)
