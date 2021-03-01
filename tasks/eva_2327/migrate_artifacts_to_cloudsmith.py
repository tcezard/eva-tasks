import argparse
import glob
import os
import re
import subprocess
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.logger import logging_config

logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)


def migrate_artifacts(python_path, cloudsmith_path, artifact_source_dir):
    # only consider directories with actual artifacts in them i.e., directories with version number names
    artifact_dirname_pattern = re.compile('[0-9]+\.[0-9]+.*')
    for dir_path, _, file_names in os.walk(artifact_source_dir):
        if artifact_dirname_pattern.match(os.path.basename(dir_path)):
            # Snapshot JARs and POMs are named in a sorted fashion but we only need the latest snapshot
            jar_file_list, pom_file_list = glob.glob(dir_path + "/*.jar"), glob.glob(dir_path + "/*.pom")
            if len(pom_file_list) > 0:
                pom_file_to_upload = sorted(pom_file_list)[-1]
                # If JAR is available for an artifact upload JAR with POM as reference
                # (ex: component libraries like accession-commons-mongodb)
                if len(jar_file_list) > 0:
                    jar_file_to_upload = sorted(jar_file_list)[-1]
                    try:
                        run_command_with_output("Migrating files {0} and {1}...".format(jar_file_to_upload,
                                                                                        pom_file_to_upload),
                                                "{0} {1} push maven ebivariation/packages {2} --pom-file={3}".format(
                                                    python_path, cloudsmith_path, jar_file_to_upload, pom_file_to_upload))
                    except subprocess.CalledProcessError as ex:
                        logger.error(ex)
                # If only POM is available, upload just the POM file
                # (ex: top-level libraries like accession-commons)
                else:
                    try:
                        run_command_with_output(
                            "Migrating files {0}...".format(pom_file_to_upload),
                            "{0} {1} push maven ebivariation/packages {2}".format(python_path, cloudsmith_path,
                                                                                  pom_file_to_upload))
                    except subprocess.CalledProcessError as ex:
                        logger.error(ex)


def main():
    description = (
        "This script migrates artifacts from JFrog artifactory to CloudSmith")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-p", "--python-path", help="Path to Python interpreter", required=True)
    parser.add_argument("-c", "--cloudsmith-path", help="Path to Cloudsmith CLI binary", required=True)
    parser.add_argument("-s", "--artifact-source-dir", help="Source directory containing the artifacts", required=True)
    args = parser.parse_args()
    migrate_artifacts(args.python_path, args.cloudsmith_path, args.artifact_source_dir)


if __name__ == "__main__":
    main()
