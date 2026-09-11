#!/usr/bin/env python
import argparse
import csv
import sys
import re

import requests
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle, resolve_variant_warehouse_db_name
from ebi_eva_internal_pyutils.mongo_utils import get_mongo_connection_handle
from ebi_eva_internal_pyutils.pg_utils import get_all_results_for_query

logger = log_cfg.get_logger(__name__)

FILES_COLLECTION = 'files_2_0'
# Project accessions are expected to look like PRJEB12345, PRJNA12345, PRJDB12345, etc. This is only
# used as a sanity check before interpolating the study id into a SQL query.
STUDY_ID_PATTERN = re.compile(r'^[A-Za-z0-9_.-]+$')
FIELD_SEPARATOR = ';'

# Public catalog of taxonomy/assembly codes for species officially released on the EVA site.
SPECIES_LIST_URL = 'https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/species/list'

REPORT_HEADER = [
    'study_id', 'inspected_db', 'assembly_accessions', 'taxonomy_ids', 'resolved_dbs',
    'db_mismatch', 'found_elsewhere', 'doc_counts_in_resolved_dbs', 'incomplete_evapro_metadata'
]

STATUS_CORRECT_NAME = 'CORRECT_NAME'
STATUS_INCORRECT_NAME_NO_STUDY = 'INCORRECT_NAME_NO_STUDY'
STATUS_INCORRECT_NAME_STUDIES_NOT_ELSEWHERE = 'INCORRECT_NAME_STUDIES_NOT_ELSEWHERE'
STATUS_INCORRECT_NAME_STUDIES_ELSEWHERE = 'INCORRECT_NAME_STUDIES_ELSEWHERE'
STATUS_DESCRIPTIONS = {
    STATUS_CORRECT_NAME: 'has the expected name: all its studies belong there',
    STATUS_INCORRECT_NAME_NO_STUDY: 'has an incorrect name and contains no study',
    STATUS_INCORRECT_NAME_STUDIES_NOT_ELSEWHERE: 'has an incorrect name and none of its misplaced '
                                                 'studies exist yet in their correct database',
    STATUS_INCORRECT_NAME_STUDIES_ELSEWHERE: 'has an incorrect name and at least one of its misplaced '
                                             'studies already exists in its correct database',
}
DATABASE_STATUS_HEADER = [
    'database', 'status', 'nb_studies', 'nb_mismatched', 'nb_found_elsewhere', 'nb_unresolved',
    'parsed_taxonomy_code', 'parsed_assembly_code', 'resolvable_via_species_list',
    'species_list_assembly_accession', 'species_list_scientific_name'
]


def get_study_ids_from_db(mongo_conn, database_name):
    """Return the sorted set of distinct study ids ("sid") found in <database_name>.files_2_0."""
    collection = mongo_conn[database_name][FILES_COLLECTION]
    study_ids = sorted(str(sid) for sid in collection.distinct('sid'))
    logger.info(f'Found {len(study_ids)} distinct study id(s) in {database_name}.{FILES_COLLECTION}: {study_ids}')
    return study_ids


def get_analyses_for_study(metadata_connection_handle, study_id):
    """
    Query EVAPRO for the distinct (analysis_accession, assembly_accession, taxonomy_id) triples
    associated with a study (project) accession.
    """
    if not STUDY_ID_PATTERN.match(study_id):
        logger.warning(f'Study id {study_id!r} does not look like a valid accession, skipping EVAPRO lookup')
        return []
    query = (
        "select distinct pa.analysis_accession, a.vcf_reference_accession, asm.taxonomy_id "
        "from project_analysis pa "
        "join analysis a on pa.analysis_accession = a.analysis_accession "
        "left join assembly_set asm on asm.assembly_set_id = a.assembly_set_id "
        f"where pa.project_accession = '{study_id}' "
        "order by pa.analysis_accession"
    )
    return get_all_results_for_query(metadata_connection_handle, query)


def study_exists_in_db(mongo_conn, database_name, study_id):
    """Return the number of documents in <database_name>.files_2_0 with the given sid."""
    collection = mongo_conn[database_name][FILES_COLLECTION]
    return collection.count_documents({'sid': study_id})


def resolve_dbs_for_analyses(metadata_connection_handle, analyses, ncbi_api_key=None):
    """
    Given a list of (analysis_accession, assembly_accession, taxonomy_id) triples, resolve the
    variant warehouse database name for each distinct (assembly_accession, taxonomy_id) pair.
    Returns a sorted list of the distinct resolved database names (assembly/taxonomy pairs missing
    either value, or that fail to resolve, are skipped).
    """
    resolved_db_per_assembly = {}
    for _, assembly_accession, taxonomy_id in analyses:
        if not assembly_accession or not taxonomy_id:
            continue
        if (assembly_accession, taxonomy_id) not in resolved_db_per_assembly:
            resolved_db_per_assembly[(assembly_accession, taxonomy_id)] = resolve_variant_warehouse_db_name(
                metadata_connection_handle, assembly_accession, taxonomy_id, ncbi_api_key=ncbi_api_key
            )
    return sorted({db for db in resolved_db_per_assembly.values() if db})


def inspect_study(mongo_conn, metadata_connection_handle, database_name, study_id, ncbi_api_key=None):
    """Build a single summary row describing where study_id (found in database_name) should live."""
    row = dict.fromkeys(REPORT_HEADER)
    row.update(study_id=study_id, inspected_db=database_name)

    analyses = get_analyses_for_study(metadata_connection_handle, study_id)
    if not analyses:
        logger.warning(f'Study {study_id} (found in {database_name}) has no analysis registered in EVAPRO')
        return row

    row['incomplete_evapro_metadata'] = any(not assembly or not taxonomy for _, assembly, taxonomy in analyses)
    row['assembly_accessions'] = FIELD_SEPARATOR.join(sorted({a for _, a, t in analyses if a}))
    row['taxonomy_ids'] = FIELD_SEPARATOR.join(sorted({str(t) for _, a, t in analyses if t}))

    resolved_dbs = resolve_dbs_for_analyses(metadata_connection_handle, analyses, ncbi_api_key=ncbi_api_key)
    if not resolved_dbs:
        logger.warning(f'Could not resolve a database name for study {study_id} (found in {database_name})')
        return row
    row['resolved_dbs'] = FIELD_SEPARATOR.join(resolved_dbs)
    row['db_mismatch'] = resolved_dbs != [database_name]

    doc_counts = {db: study_exists_in_db(mongo_conn, db, study_id) for db in resolved_dbs}
    row['doc_counts_in_resolved_dbs'] = FIELD_SEPARATOR.join(f'{db}:{count}' for db, count in doc_counts.items())

    if row['db_mismatch']:
        row['found_elsewhere'] = any(count > 0 for db, count in doc_counts.items() if db != database_name)
        if row['found_elsewhere']:
            logger.warning(
                f'Study {study_id} found in both {database_name} (should not be there) and its correct '
                f'database(s): {row["doc_counts_in_resolved_dbs"]}'
            )
        else:
            logger.warning(
                f'Study {study_id} found in {database_name} but should be in {row["resolved_dbs"]}, '
                f'which does not currently contain it'
            )
    else:
        row['found_elsewhere'] = False
        logger.info(
            f'Study {study_id} is correctly placed in {database_name} ({row["doc_counts_in_resolved_dbs"]})'
        )
    return row


def inspect_database(mongo_conn, metadata_connection_handle, database_name, ncbi_api_key=None):
    """Inspect every study present in database_name and return one report row per study."""
    study_ids = get_study_ids_from_db(mongo_conn, database_name)
    return [
        inspect_study(mongo_conn, metadata_connection_handle, database_name, study_id, ncbi_api_key=ncbi_api_key)
        for study_id in study_ids
    ]


def inspect_databases(mongo_conn, metadata_connection_handle, database_names, ncbi_api_key=None):
    """Inspect each database in database_names in turn and return the concatenated report rows."""
    rows = []
    for database_name in database_names:
        rows.extend(inspect_database(mongo_conn, metadata_connection_handle, database_name,
                                     ncbi_api_key=ncbi_api_key))
    return rows


def get_species_list(species_list_url=SPECIES_LIST_URL):
    """
    Fetch the public EVA species list, used to independently verify a database name's taxonomy and
    assembly codes.
    """
    try:
        response = requests.get(species_list_url, timeout=30)
        response.raise_for_status()
        return response.json()['response'][0]['result']
    except Exception as e:
        logger.warning(f'Could not fetch the public species list from {species_list_url}: {e}')
        return None


def parse_db_name(database_name):
    """
    Split a variant warehouse database name into (taxonomy_code, assembly_code).
    """
    parts = database_name.split('_')
    if len(parts) < 3 or parts[0] != 'eva':
        return None, None
    return parts[1], '_'.join(parts[2:])


def check_db_name_against_species_list(database_name, species_list):
    """
    Cross-check a database name against the public EVA species list
    """
    taxonomy_code, assembly_code = parse_db_name(database_name)
    result = {
        'parsed_taxonomy_code': taxonomy_code, 'parsed_assembly_code': assembly_code,
        'resolvable_via_species_list': None, 'species_list_assembly_accession': None,
        'species_list_scientific_name': None
    }
    if taxonomy_code is None:
        logger.warning(
            f'{database_name} does not look like "eva_<taxonomy_code>_<assembly_code>", cannot check '
            f'it against the public species list'
        )
        return result
    if species_list is None:
        return result

    matches = [r for r in species_list
              if r.get('taxonomyCode') == taxonomy_code and r.get('assemblyCode') == assembly_code]
    result['resolvable_via_species_list'] = bool(matches)
    if matches:
        # Keep the most recent patch when several assemblies share the same code
        best_match = max(matches, key=lambda r: r.get('assemblyAccession', ''))
        result['species_list_assembly_accession'] = best_match.get('assemblyAccession')
        result['species_list_scientific_name'] = best_match.get('taxonomyScientificName')
    return result


def classify_database(database_name, rows_for_db, species_list=None):
    """
    Determine the overall status of one inspected database from its per-study report rows.
    The database name is also checked against the public species list via API
    """
    mismatched = [r for r in rows_for_db if r['db_mismatch']]
    unresolved = [r for r in rows_for_db if r['db_mismatch'] is None]
    found_elsewhere = [r for r in mismatched if r['found_elsewhere']]

    if not rows_for_db:
        status = STATUS_INCORRECT_NAME_NO_STUDY
    elif not mismatched:
        status = STATUS_CORRECT_NAME
    elif found_elsewhere:
        status = STATUS_INCORRECT_NAME_STUDIES_ELSEWHERE
    else:
        status = STATUS_INCORRECT_NAME_STUDIES_NOT_ELSEWHERE

    if unresolved:
        logger.warning(
            f'{database_name}: {len(unresolved)} stud(ies) could not be checked against EVAPRO and were '
            f'left out of the status determination'
        )

    name_check = check_db_name_against_species_list(database_name, species_list)
    if name_check['resolvable_via_species_list'] is False and status == STATUS_CORRECT_NAME:
        logger.warning(
            f'{database_name}: EVAPRO-based check found no mismatch (status={status}) but the public '
            f'species list has no entry for taxonomy code "{name_check["parsed_taxonomy_code"]}" and '
            f'assembly code "{name_check["parsed_assembly_code"]}" - the two checks disagree, please '
            f'investigate by hand'
        )
    elif name_check['resolvable_via_species_list'] is True and status != STATUS_CORRECT_NAME:
        logger.warning(
            f'{database_name}: EVAPRO-based check reports {status} but the public species list does '
            f'resolve this name to {name_check["species_list_scientific_name"]} '
            f'({name_check["species_list_assembly_accession"]}) - the two checks disagree, please '
            f'investigate by hand'
        )

    return {
        'database': database_name,
        'status': status,
        'nb_studies': len(rows_for_db),
        'nb_mismatched': len(mismatched),
        'nb_found_elsewhere': len(found_elsewhere),
        'nb_unresolved': len(unresolved),
        **name_check,
    }


def summarize_database_statuses(rows, database_names, species_list=None):
    """Classify each database in database_names from the per-study rows built by inspect_databases."""
    return [
        classify_database(database_name, [r for r in rows if r['inspected_db'] == database_name],
                         species_list=species_list)
        for database_name in database_names
    ]


def write_report(rows, output):
    writer = csv.DictWriter(output, fieldnames=REPORT_HEADER, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)


def write_database_status_report(statuses, output):
    writer = csv.DictWriter(output, fieldnames=DATABASE_STATUS_HEADER, delimiter='\t')
    writer.writeheader()
    writer.writerows(statuses)


def main():
    parser = argparse.ArgumentParser(
        description='Inspect one or more variant warehouse databases for studies that might belong in a '
                    'different, correctly named database.'
    )
    parser.add_argument('--database-names', required=True, nargs='+',
                        help='Name(s) of the (suspect) database(s) in the variant warehouse to inspect')
    parser.add_argument('--settings-xml-file', required=True,
                        help='Path to the private settings XML file used to connect to mongodb and EVAPRO')
    parser.add_argument('--profile', default='production', help='Maven profile to use (default: production)')
    parser.add_argument('--eutils-api-key',
                        help='NCBI eutils API key, used to resolve assemblies not yet known to EVAPRO')
    parser.add_argument('--species-list-url', default=SPECIES_LIST_URL,
                        help='URL of the public EVA species list used to cross-check database names '
                             '(default: %(default)s)')
    parser.add_argument('--skip-species-list-check', action='store_true',
                        help='Skip the public species list cross-check (e.g. if offline)')
    parser.add_argument('--output', help='Path to the TSV report to write (default: stdout)')
    args = parser.parse_args()
    log_cfg.add_stdout_handler()

    with get_mongo_connection_handle(args.profile, args.settings_xml_file) as mongo_conn, \
            get_metadata_connection_handle(args.profile, args.settings_xml_file) as metadata_connection_handle:
        rows = inspect_databases(mongo_conn, metadata_connection_handle, args.database_names,
                                 ncbi_api_key=args.eutils_api_key)

    write_report(rows, sys.stdout)

    species_list = None if args.skip_species_list_check else get_species_list(args.species_list_url)
    statuses = summarize_database_statuses(rows, args.database_names, species_list=species_list)
    logger.info('Database status summary:')
    for db_status in statuses:
        logger.info(
            f"{db_status['database']}: {db_status['status']} "
            f"({STATUS_DESCRIPTIONS[db_status['status']]}) - studies={db_status['nb_studies']}, "
            f"mismatched={db_status['nb_mismatched']}, found_elsewhere={db_status['nb_found_elsewhere']}, "
            f"unresolved={db_status['nb_unresolved']}, "
            f"resolvable_via_species_list={db_status['resolvable_via_species_list']}"
        )
    if args.output:
        with open(args.output, 'w', newline='') as output_file:
            write_database_status_report(statuses, output_file)
        logger.info(f'Report written to {args.output}')
    else:
        write_database_status_report(statuses, sys.stdout)


if __name__ == '__main__':
    main()
