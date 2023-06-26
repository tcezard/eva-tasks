import logging
import os
import shutil
from argparse import ArgumentParser

from ebi_eva_common_pyutils import command_utils
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)
log_cfg.set_log_level(logging.INFO)


def main():
    argparse = ArgumentParser(description='COPY logs from clustering to logs directory')
    argparse.add_argument('--projects_dir_path', required=True, type=str, help='The Path to projects dir')
    argparse.add_argument('--project', required=False, type=str,
                          help='Project for which logs needs to be copied/deleted')
    argparse.add_argument('--copy_logs', action='store_true', default=False,
                          help='Copy logs from clustering directory to logs directory')
    argparse.add_argument('--delete_clustering_logs', action='store_true', default=False,
                          help='Delete logs from clustering directory')

    args = argparse.parse_args()

    if not (args.copy_logs or args.delete_clustering_logs):
        argparse.error(f'No action requested, add --copy_logs or --delete_clustering_logs')

    if args.project:
        project_list = [args.project]
    else:
        project_list = [f for f in os.listdir(args.projects_dir_path) if f[:5] == 'PRJEB']

    for project in project_list:
        logger.info(f'Processing project {project}')
        project_clustering_logs_dir = os.path.join(args.projects_dir_path, project, '53_clustering', 'logs')

        if os.path.exists(project_clustering_logs_dir):
            if args.copy_logs:
                project_logs_dir = os.path.join(args.projects_dir_path, project, '00_logs')
                if os.listdir(project_clustering_logs_dir):
                    logger.info(f'Copying Files for project {project}')

                    copy_command = f'cp {project_clustering_logs_dir}/* {project_logs_dir}'
                    command_utils.run_command_with_output('copying log files', copy_command)
                else:
                    logger.warning(f'Logs directory is empty for project {project}')

            if args.delete_clustering_logs:
                shutil.rmtree(project_clustering_logs_dir)
        else:
            logger.warning(f'Clustering logs directory does not exist for project {project}')


if __name__ == "__main__":
    main()
