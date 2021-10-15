from argparse import ArgumentParser
from datetime import datetime

import psycopg2
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query


def serialise_tuple(t):
    res = []
    for e in t:
        if e is None:
            res.append("null")
        elif isinstance(e, str) or isinstance(e, datetime):
            res.append(f"'{e}'")
        else:
            res.append(f"{e}")
    return '(' + ', '.join(res) + ')'


class RestoreProject:

    def __init__(self, prod_conn, dev_conn):
        self.prod_conn = prod_conn
        self.dev_conn = dev_conn
        self.dev_cursor = self.dev_conn.cursor()

    def restore_table(self, table_name, where_clause):
        query = f"select * from {table_name} where {where_clause}"
        rows = get_all_results_for_query(self.prod_conn, query)
        for row in rows:
            query = f'insert into {table_name} VALUES ' + serialise_tuple(row)
            print(query)
            self.dev_cursor.execute(query)

    def restore(self, project_accession):
        self.restore_table('project', where_clause=f"project_accession='{project_accession}'")

        query = "select dbxref_id from project_dbxref where project_accession='{project_accession}'"
        rows = get_all_results_for_query(self.prod_conn, query)
        for dbxref_id, in rows:
            self.restore_table('dbxref', where_clause=f'where dbxref_id={dbxref_id}')
        self.restore_table('project_dbxref', where_clause=f"project_accession='{project_accession}'")
        self.restore_table('project_taxonomy', where_clause=f"project_accession='{project_accession}'")

        query = f"select submission_id from project_ena_submission where project_accession='{project_accession}'"
        rows = get_all_results_for_query(self.prod_conn, query)
        for submission_id, in rows:
            self.restore_table('submission', where_clause=f"submission_id={submission_id}")

        query = f"select analysis_accession from project_analysis where project_accession='{project_accession}'"
        rows = get_all_results_for_query(self.prod_conn, query)
        for analysis_accession, in rows:
            query = f"select file_id from file join analysis_file using(file_id) where analysis_accession='{analysis_accession}'"
            rows = get_all_results_for_query(self.prod_conn, query)
            for file_id, in rows:
                self.restore_table('file', where_clause=f"file_id={file_id}")
                self.restore_table('browsable_file', where_clause=f"file_id={file_id}")
            self.restore_table('analysis', where_clause=f"analysis_accession='{analysis_accession}'")
            self.restore_table('project_analysis', where_clause=f"analysis_accession='{analysis_accession}'")
            self.restore_table('analysis_file', where_clause=f"analysis_accession='{analysis_accession}'")
            self.restore_table('analysis_sequence', where_clause=f"analysis_accession='{analysis_accession}'")
            self.restore_table('analysis_submission', where_clause=f"analysis_accession='{analysis_accession}'")
            self.restore_table('analysis_experiment_type', where_clause=f"analysis_accession='{analysis_accession}'")
            self.restore_table('analysis_platform', where_clause=f"analysis_accession='{analysis_accession}'")

        self.restore_table('project_ena_submission', where_clause=f"project_accession='{project_accession}'")

        query = f"select eload_id from project_eva_submission where project_accession='{analysis_accession}'"
        rows = get_all_results_for_query(self.prod_conn, query)
        for eload_id, in rows:
            self.restore_table('eva_submission', where_clause=f"eva_submission_id={eload_id}")

        self.dev_conn.commit()


def main():
    parser = ArgumentParser()
    parser.add_argument('--project_accession', required=True)
    parser.add_argument('--production_uri', required=True)
    parser.add_argument('--production_user', required=True)
    parser.add_argument('--production_pass', required=True)
    parser.add_argument('--development_uri', required=True)
    parser.add_argument('--development_user', required=True)
    parser.add_argument('--development_pass', required=True)
    args = parser.parse_args()

    dev_conn = psycopg2.connect(args.development_uri, user=args.development_user, password=args.development_pass)
    prod_conn = psycopg2.connect(args.production_uri, user=args.production_user, password=args.production_pass)

    RestoreProject(prod_conn, dev_conn).restore(args.project_accession)

    prod_conn.close()
    dev_conn.close()


if __name__ == '__main__':
    main()

