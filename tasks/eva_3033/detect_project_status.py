import argparse
import glob
import gzip
import os.path

from cached_property import cached_property

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query, execute_query


logger = logging_config.get_logger(__name__)
logging_config.add_stderr_handler()


def detect_all_public_project(maven_config, maven_profile):
    with get_metadata_connection_handle(maven_profile, maven_config) as pg_conn:
        query = f"select project_accession from evapro.project where ena_status=4"
        public_projects = [project for project, in get_all_results_for_query(pg_conn, query)]

    return public_projects


def process_projects(projects, maven_config, maven_profile, noah_project_dir, codon_project_dir):
    for project in projects:
        detector = ProjectStatusDetector(project, maven_config, maven_profile, noah_project_dir, codon_project_dir)
        detector.detect_project_status()


class ProjectStatusDetector:

    def __init__(self, project, maven_config, maven_profile, noah_project_dir, codon_project_dir):
        self.project = project
        self.maven_config = maven_config
        self.maven_profile = maven_profile
        self.noah_project_dir = noah_project_dir
        self.codon_project_dir = codon_project_dir

    @cached_property
    def metadata_conn(self):
        return get_metadata_connection_handle(self.maven_profile, self.maven_config)

    @cached_property
    def mongo_conn(self):
        return get_mongo_connection_handle(self.maven_profile, self.maven_config)

    def detect_project_status(self):
        for analysis, source_assembly, taxonomy, filenames in self.project_information():
            accessioning_status = remapping_status = clustering_status = target_assembly = 'Not found'
            if taxonomy and taxonomy != 9606:
                list_ssid = self.check_accessioning_was_done(analysis, filenames)
                accessioning_status = 'Done' if len(list_ssid) > 0 else 'Pending'
                target_assembly = self.find_current_target_assembly_for(taxonomy)
                remapping_status = 'Required' if source_assembly != target_assembly else 'Not_required'
                if source_assembly != target_assembly:
                    if self.check_remapping_was_done(source_assembly, target_assembly, list_ssid):
                        remapping_status = 'Done'
                    assembly = target_assembly
                else:
                    assembly = source_assembly

                clustering_status = 'Done' if self.check_clustering_was_done(assembly, list_ssid) else 'Pending'
            else:
                if not taxonomy:
                    logger.error( f'Project {self.project}:{analysis} has no taxonomy associated and the metadata '
                              f'should be checked.')
            print('\t'.join([self.project, str(analysis), str(taxonomy), str(source_assembly), str(target_assembly),
                             accessioning_status, remapping_status, clustering_status]))

    def project_information(self):
        """Retrieve project information from the metadata. Information retrieve include
        the analysis and associated taxonomy, genome and file names that are included in this project."""
        query = (
            "select pa.project_accession, pa.analysis_accession, a.vcf_reference_accession, at.taxonomy_id, f.filename "
            "from project_analysis pa "
            "join analysis a on pa.analysis_accession=a.analysis_accession "
            "join assembly_set at on at.assembly_set_id=a.assembly_set_id "
            "join analysis_file af on af.analysis_accession=a.analysis_accession "
            "join file f on f.file_id=af.file_id "
            f"where f.file_type='VCF' and pa.project_accession='{self.project}'"
            "order by pa.project_accession, pa.analysis_accession")
        filenames = []
        current_analysis = current_assembly = current_tax_id = None
        for project, analysis, assembly, tax_id, filename in get_all_results_for_query(self.metadata_conn, query):
            if analysis != current_analysis:
                if current_analysis:
                    yield current_analysis, current_assembly, current_tax_id, filenames
                current_analysis = analysis
                current_assembly = assembly
                current_tax_id = tax_id
                filenames = []
            filenames.append(filename)
        yield current_analysis, current_assembly, current_tax_id, filenames

    def check_accessioning_was_done(self, analysis, filenames):
        """
        Check that an accessioning file can be found in either noah or codon (assume access to both filesystem)
        It parses and provide a 1000 submitted variant accessions from that project.
        """
        accessioning_reports = self.get_accession_reports_for_study()
        accessioned_filenames = [self.get_accession_file(f) for f in filenames]
        if not accessioning_reports:
            return 0, []
        if len(accessioning_reports) == 1:
            accessioning_report = accessioning_reports[0]
        elif len([r for r in accessioning_reports if r in accessioned_filenames]) == 1:
            accessioning_report = [r for r in accessioning_reports if os.path.basename(r) in accessioned_filenames][0]
        elif len([r for r in accessioning_reports if analysis in r]) == 1:
            accessioning_report = [r for r in accessioning_reports if analysis in r][0]
        elif accessioning_reports:
            logger.warning(f'Assume all accessioning report are from project {self.project}:{analysis} and only use the first one.')
            accessioning_report = accessioning_reports[0]
        else:
            logger.error(f'Cannot assign accessioning report to project {self.project} analysis {analysis} for files {accessioning_reports}')
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
        Read the accessioning report to retrieve the first 1000 Submitted variant accessions
        """
        no_of_ss_ids_in_file = 0
        first_1000_ids = []
        with gzip.open(path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                elif line.split("\t")[2].startswith("ss"):
                    no_of_ss_ids_in_file = no_of_ss_ids_in_file + 1
                    if no_of_ss_ids_in_file <= 1000:
                        first_1000_ids.append(int(line.split("\t")[2][2:]))
                    else:
                        break
        return first_1000_ids

    def get_accession_reports_for_study(self):
        """
        Given a study, find the accessioning report path for that study on both noah and codon
        look for a file ending with accessioned.vcf.gz (assuming only one file with the given pattern will be present)
        if more than one files are present, take the first one
        """

        noah_files = glob.glob(os.path.join(self.noah_project_dir, self.project, '60_eva_public', '*accessioned.vcf.gz'))
        codon_files = glob.glob(os.path.join(self.codon_project_dir, self.project, '60_eva_public', '*accessioned.vcf.gz'))
        accessioning_reports = noah_files + codon_files
        if not accessioning_reports:
            logger.error(f"Could not find any file in Noah or Codon for Study {self.project}")

        return accessioning_reports

    def check_remapping_was_done(self, source_assembly, target_assembly, list_ssid):
        # query the tracker to filter out things that were never remapped
        query = (f'select release_version, remapping_status, remapping_end ' 
                 f'from eva_progress_tracker.remapping_tracker '
                 f"where origin_assembly_accession='{source_assembly}' and assembly_accession='{target_assembly}'")
        remapped_candidate = None
        for release_version, remapping_status, remapping_end in get_all_results_for_query(self.metadata_conn, query):
            remapped_candidate = True
        if remapped_candidate:
            ss_variants = self.find_submitted_variant_in_assembly(target_assembly, list_ssid)
            logger.info(f'Found {len(ss_variants)} variants out of {len(list_ssid)} in {target_assembly}')
            return len(ss_variants) > 0
        return False

    def check_clustering_was_done(self, assembly, list_ssid):
        ss_variants = self.find_submitted_variant_in_assembly(assembly, list_ssid)
        logger.info(f'Found {len(ss_variants)} variants out of {len(list_ssid)} in {assembly}')
        return len([ss_variant for ss_variant in ss_variants if 'rs' in ss_variant]) > 0

    def find_submitted_variant_in_assembly(self, assembly, list_ssid):
        filters = {'seq': assembly, 'accession': {'$in': list_ssid}}
        cursor = self.mongo_conn['eva_accession_sharded']['submittedVariantEntity'].find(filters)
        variants = []
        for variant in cursor:
            variants.append(variant)
        return variants

    def find_current_target_assembly_for(self, taxonomy):
        query = f"select assembly_id from evapro.supported_assembly_tracker where taxonomy_id={taxonomy} and current=true"
        assemblies = []
        for asm, in get_all_results_for_query(self.metadata_conn, query):
            assemblies.append(asm)
        assert len(assemblies) < 2, f'Multiple target assemblies found for taxonomy {taxonomy}'
        if assemblies:
            return assemblies[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Provide the processing status of the EVA projects')
    parser.add_argument("--projects",  nargs='*', default=None,
                        help="project you was the status of. If not specified then all projects are checked",)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)

    parser.add_argument("--noah-prj-dir", help="path to the project directory in noah", required=True)
    parser.add_argument("--codon-prj-dir", help="path to the project directory in codon", required=True)
    parser.add_argument("--profile", choices=('localhost', 'development', 'production_processing'),
                        help="Profile to decide which environment should be used for making entries", required=True)

    args = parser.parse_args()
    if args.projects:
        projects = args.projects
    else:
        projects = detect_all_public_project(args.private_config_xml_file, args.profile)
    process_projects(projects, args.private_config_xml_file, args.profile, args.noah_prj_dir, args.codon_prj_dir)
