import argparse
import glob
import gzip
import os.path
from collections import defaultdict
from urllib.parse import urlsplit

import cached_property
import psycopg2
from ebi_eva_common_pyutils.config_utils import get_accession_pg_creds_for_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query, execute_query
from ebi_eva_common_pyutils.taxonomy.taxonomy import get_scientific_name_from_ensembl

from genome_target_tracker import get_tax_latest_asm_from_eva, add_assembly_to_accessioned_assemblies

logger = logging_config.get_logger(__name__)


class ClusteringDetector:

    def __init__(self, maven_config, maven_profile, noah_project_dir, codon_project_dir):
        self.maven_config = maven_config
        self.maven_profile = maven_profile
        self.noah_project_dir = noah_project_dir
        self.codon_project_dir = codon_project_dir


    @cached_property
    def metadata_conn(self):
        return get_metadata_connection_handle(self.maven_profile, self.maven_profile)

    @cached_property
    def mongo_conn(self):
        return get_mongo_connection_handle(self.maven_profile, self.maven_profile)

    def detect_project_needing_clustering(self):
        # From project get source assembly and taxonomy
        # From taxonomy find target assembly
        # if target assembly != source assembly then it should have been remapped
        #   Check that it has been remapped
        # source assembly = target assembly
        # Check that the submitted variants have been clustered

        for project, analysis, source_assembly, taxonomy, filenames in self.get_project_information():
            target_assembly = self.find_current_target_assembly_for(taxonomy)
            nb_ss_id, list_ssid = self.check_accessioning_was_done(project, analysis, filenames)
            if nb_ss_id == 0:
                logger.error(f'Project {project}:{analysis} was not accessionned')
            if source_assembly != target_assembly:
                logger.info(f'Remapping required for project {project}:{analysis} from {source_assembly} to {target_assembly}')
                self.check_remapping_was_done(source_assembly, target_assembly, list_ssid)
                assembly = target_assembly
            else:
                assembly = source_assembly
            self.check_clustering_was_done(assembly, list_ssid)

    def check_accessioning_was_done(self, project, analysis, filenames):
        """
        For each taxonomy and assembly, get the total number of ssids to remap, by reading the accessioning report
        for all the studies which needs to be remapped
        """
        accessioning_reports = self.get_accession_reports_for_study(project)
        accessioned_filenames = [self.get_accession_file(f) for f in filenames]
        if not accessioning_reports:
            return 0, []
        if len(accessioning_reports) == 1:
            accessioning_report = accessioning_reports[0]
        elif len([r for r in accessioning_reports if r in accessioned_filenames]) == 1:
            accessioning_report = [r for r in accessioning_reports if os.path.basename(r) in accessioned_filenames][0]
        elif len([r for r in accessioning_reports if analysis in r]) == 1:
            accessioning_report = [r for r in accessioning_reports if analysis in r][0]
        else:
            logger.error(f'Cannot assign accessioning report to project {project} analysis {analysis} for files {accessioning_reports}')
            accessioning_report = None
        if not accessioning_report:
            return 0, []
        return self.get_accessioning_info_from_file(accessioning_report)

    def get_accession_file(self, filename):
        basefile = ''
        if filename.endswith('.vcf.gz'):
            basefile = filename[:-7]
        elif filename.endswith('.vcf'):
            basefile = filename[:-4]
        return basefile + '.accessioned.vcf.gz '

    def get_accessioning_info_from_file(self, path):
        """
        Read a file to figure out the number of ss ids in that file
        skip line start with #
        read line where 3rd entry starts with ss (assuming all files have same sequence of headers and ss id is always 3rd)
        """
        no_of_ss_ids_in_file = 0
        first_1000_ids = []
        with gzip.open(path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                elif line.split("\t")[2].startswith("ss"):
                    no_of_ss_ids_in_file = no_of_ss_ids_in_file + 1
                    if no_of_ss_ids_in_file < 1000:
                        first_1000_ids.append(line.split("\t")[2])

        logger.info(f"No of ss ids in file {path} : {no_of_ss_ids_in_file}")
        return no_of_ss_ids_in_file, first_1000_ids


    def get_accession_reports_for_study(self, project):
        """
        Given a study, find the accessioning report path for that study on both noah and codon
        look for a file ending with accessioned.vcf.gz (assuming only one file with the given pattern will be present)
        if more than one files are present, take the first one
        """

        noah_files = glob.glob(f"{self.noah_project_dir}/{project}/60_eva_public/*accessioned.vcf.gz")
        codon_files = glob.glob(f"{self.noah_project_dir}/{project}/60_eva_public/*accessioned.vcf.gz")
        accessioning_reports = noah_files + codon_files
        if not accessioning_reports:
            logger.error(f"Could not find any file in Noah or Codon for Study {project}")

        return accessioning_reports

    def check_remapping_was_done(self, source_assembly, target_assembly, list_ssid):
        # query the tracker to filter out things that were never remapped
        query = (f'select release_version, remapping_status, remapping_end ' 
                 f'from eva_progress_tracker.remapping_tracker;'
                 f"where origin_assembly_accession='{source_assembly}' and assembly_accession='{target_assembly}'")
        remapped_candidate = None
        for release_version, remapping_status, remapping_end in get_all_results_for_query(self.metadata_conn, query):
            remapped_candidate = True
        if remapped_candidate:
            ss_variants = self.find_submitted_variant_in_assembly(target_assembly, list_ssid)
            logger.info(f'Found {len(ss_variants)} variants out of {len(list_ssid)} in {target_assembly}')

    def check_clustering_was_done(self, assembly, list_ssid):
        ss_variants = self.find_submitted_variant_in_assembly(assembly, list_ssid)
        logger.info(f'Found {len(ss_variants)} variants out of {len(list_ssid)} in {assembly}')

    def find_submitted_variant_in_assembly(self, assembly, list_ssid):
        list_accession = [int(ss.strip('s')) for ss in list_ssid]
        cursor = self.mongo_conn['eva_accessioned_sharded']['submittedVariantEntity'].find({'seq': assembly, 'accession': {'$in': list_accession}})
        variants = []
        for variant in cursor:
            variants.append(variant)
        return variants

    def get_project_information(self):
        query = ("select pa.project_accession, pa.analysis_accession, a.vcf_reference_accession, at.taxonomy_id, f.filename "
                 "from project_analysis pa "
                 "join analysis a on pa.analysis_accession=a.analysis_accession "
                 "join assembly_set at on at.assembly_set_id=a.assembly_set_id "
                 "join analysis_file af on af.analysis_accession=a.analysis_accession "
                 "join file f on f.file_id=af.file_id "
                 "where f.file_type='VCF'"
                 "order by pa.project_accession, pa.analysis_accession")
        filenames = []
        current_project = current_analysis = current_assembly = curren_tax_id = None
        for project, analysis, assembly, tax_id, filename in get_all_results_for_query(self.metadata_conn, query):
            if analysis != current_analysis:
                if analysis:
                    yield current_project, current_analysis, current_assembly, curren_tax_id, filenames
                current_project = project
                current_analysis = analysis
                current_assembly = assembly
                curren_tax_id = tax_id
                filenames = []
            filenames.append(filename)
        yield current_project, current_analysis, current_assembly, curren_tax_id, filenames


    def find_current_target_assembly_for(self, taxonomy):
        query = f"select assembly_id from evapro.supported_assembly_tracker where taxonomy_id={taxonomy} and current=true"
        assemblies = []
        for asm in get_all_results_for_query(self.metadata_conn, query):
            assemblies.append(asm)
        assert len(assemblies) < 2, f'Multiple target assemblies found for taxonomy {taxonomy}'
        if assemblies:
            return assemblies[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare remapping table for release', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--noah-prj-dir", help="path to the project directory in noah", required=True)
    parser.add_argument("--codon-prj-dir", help="path to the project directory in codon", required=True)
    parser.add_argument("--profile", choices=('localhost', 'development', 'production_processing'),
                        help="Profile to decide which environment should be used for making entries", required=True)

    args = parser.parse_args()
    detector = ClusteringDetector(args.private_config_xml_file, args.profile, args.noah_prj_dir, args.codon_prj_dir)
    detector.detect_project_needing_clustering()
