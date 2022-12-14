import argparse
import datetime
import json
import os.path
import atexit
from collections import defaultdict

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.ncbi_utils import get_ncbi_assembly_dicts_from_term, get_ncbi_taxonomy_dicts_from_ids
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query

cached_dict = {}
cached_file = 'cached_data.json'
if os.path.exists(cached_file):
    with open(cached_file) as open_file:
        cached_dict = json.load(open_file)
if 'assembly_dicts' not in cached_dict: cached_dict['assembly_dicts'] = {}
if 'taxonomy_dicts' not in cached_dict: cached_dict['taxonomy_dicts'] = {}

def exit_handler():
    with open(cached_file, 'w') as open_file:
        json.dump(cached_dict, open_file)


atexit.register(exit_handler)



logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)

def cached_get_key_from(assembly, key):
    if assembly not in cached_dict['assembly_dicts']:
        cached_dict['assembly_dicts'][assembly] = get_ncbi_assembly_dicts_from_term(assembly)
    assembly_dicts = cached_dict['assembly_dicts'][assembly]
    values = set([d.get(key) for d in assembly_dicts])
    if len(values) > 1:
        # Only keep the one that have the assembly accession as a synonymous and check again
        values = set([d.get('key') for d in assembly_dicts
                      if assembly in d['synonym'].values() or assembly == d['assemblyaccession']])
    if len(values) != 1:
        raise ValueError(f"Cannot resolve assembly's {key} for assembly {assembly} in NCBI. "
                         f'Found {",".join([str(a) for a in values])}')
    return values.pop()


def cached_get_taxonomy_from(assembly):
    return int(cached_get_key_from(assembly, 'taxid'))


def cached_get_assembly_name_from(assembly):
    return cached_get_key_from(assembly, 'assemblyname')


def cached_get_scientific_name(taxid):
    if taxid not in cached_dict['taxonomy_dicts']:
        cached_dict['taxonomy_dicts'][taxid] = get_ncbi_taxonomy_dicts_from_ids([str(taxid)])
    taxonomy_dicts = cached_dict['taxonomy_dicts'][taxid]
    scientific_names = set([d.get('scientificname') for d in taxonomy_dicts])
    if len(scientific_names) != 1:
        raise ValueError(f"Cannot resolve taxonomy's name for taxonomy_id {taxid} in NCBI. "
                         f'Found {",".join([str(a) for a in scientific_names])}')
    return scientific_names.pop()


def detect_species_with_cross_species_target(maven_config, maven_profile):
    with get_metadata_connection_handle(maven_profile, maven_config) as pg_conn:
        query = f"select taxonomy_id,assembly_id, current from evapro.supported_assembly_tracker;"
        taxonomy_to_assemblies = defaultdict(list)
        taxonomy_to_current_assembly = {}

        for taxonomy, assembly, current in get_all_results_for_query(pg_conn, query):
            taxonomy_to_assemblies[taxonomy].append(assembly)
            if current:
                taxonomy_to_current_assembly[taxonomy] = assembly
        clustered_assemblies_to_dates = defaultdict(list)

        # Get the assemblies where clustering for whole species was performed
        query = 'select assembly_accession, clustering_start from eva_progress_tracker.clustering_release_tracker where should_be_clustered=true;'
        for assembly, clustering_start in get_all_results_for_query(pg_conn, query):
            clustered_assemblies_to_dates[assembly].append(clustering_start)

        # Get the assemblies where clustering for new studies was performed
        query = 'select assembly_accession, remapping_start from eva_progress_tracker.remapping_tracker where study_accessions is not null;'
        for assembly, remapping_start in get_all_results_for_query(pg_conn, query):
            clustered_assemblies_to_dates[assembly].append(remapping_start)

    for taxonomy in taxonomy_to_current_assembly:
        for assembly in taxonomy_to_assemblies[taxonomy]:
            taxonomy_from_assembly = cached_get_taxonomy_from(assembly)
            # No difference between the taxonomy and the assembly
            if taxonomy == taxonomy_from_assembly:
                continue
            # taxonomy are different but the target assembly is the same as the current one for this
            if assembly == taxonomy_to_current_assembly.get(taxonomy_from_assembly):
                # logger.info(f'For {taxonomy} target is {taxonomy_to_assembly[taxonomy]} which comes from {taxonomy_from_assembly} with same target {taxonomy_to_assembly.get(taxonomy_from_assembly)}')
                continue
            # taxonomy are different and assembly is different from the current target
            # but the target's taxonomy has itself no remapping target
            if taxonomy_to_current_assembly.get(taxonomy_from_assembly) is None:
                # logger.info(f'For {taxonomy} target is {taxonomy_to_assembly[taxonomy]} which comes from {taxonomy_from_assembly} with same target {taxonomy_to_assembly.get(taxonomy_from_assembly)}')
                continue

            # Check that if a clustering was ever performed on this assembly
            if assembly not in clustered_assemblies_to_dates:
                continue

            species_name = cached_get_scientific_name(taxonomy)
            species_name_from_assembly = cached_get_scientific_name(taxonomy_from_assembly)
            all_dates = ','.join([d.date().isoformat() if d else 'Unknown' for d in clustered_assemblies_to_dates[assembly]])
            logger.warning(f'For {taxonomy} ({species_name}) and target {assembly} '
                           f'({cached_get_assembly_name_from(assembly)}), the assembly comes from '
                           f'{taxonomy_from_assembly} ({species_name_from_assembly}) and its current is target '
                           f'{taxonomy_to_current_assembly.get(taxonomy_from_assembly)} '
                           f'({cached_get_assembly_name_from(taxonomy_to_current_assembly.get(taxonomy_from_assembly))}). '
                           f'Clustering happened on {all_dates}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Search for species where the target assembly is from a different taxonomy')
    parser.add_argument("--maven_config", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--maven_profile", choices=('localhost', 'development', 'production_processing'),
                        help="Profile to decide which environment should be used for making entries", required=True)

    args = parser.parse_args()
    detect_species_with_cross_species_target(args.maven_config, args.maven_profile)
