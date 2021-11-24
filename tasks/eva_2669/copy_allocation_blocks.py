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


def main():
    parser = ArgumentParser()
    parser.add_argument('--production_uri', required=True)
    parser.add_argument('--production_user', required=True)
    parser.add_argument('--production_pass', required=True)
    parser.add_argument('--development_uri', required=True)
    parser.add_argument('--development_user', required=True)
    parser.add_argument('--development_pass', required=True)
    parser.add_argument('--apply_queries', default=False, action='store_true')
    args = parser.parse_args()

    dev_conn = psycopg2.connect(args.development_uri, user=args.development_user, password=args.development_pass)
    prod_conn = psycopg2.connect(args.production_uri, user=args.production_user, password=args.production_pass)
    prod_cursor = prod_conn.cursor()

    # Get the newest existing block from prod
    query = """select max("first_value") from public.contiguous_id_blocks where category_id ='rs' """
    largest_first_value_in_prod, = get_all_results_for_query(prod_conn, query)[0]

    # Fill in with incomplete blocks until reaching the first block used in dev
    for first_value in range(largest_first_value_in_prod + 100000, 3166700000, 100000):
        query = 'INSERT INTO public.contiguous_id_blocks VALUES ' + \
                f"('instance-6', 'rs', {first_value}, {first_value}, {first_value + 999999})"
        print(query)
        if args.apply_queries:
            prod_cursor.execute(query)

    # Retrieve the blocks used in dev
    query = (f'select application_instance_id, category_id, "first_value", last_committed , "last_value" '
             f'from public.contiguous_id_blocks '
             f"where category_id ='rs' "
             f"and application_instance_id = 'instance-6'"
             f"and first_value >= 3166700000;")
    rows = get_all_results_for_query(dev_conn, query)

    # Write them out in prod
    for row in rows:
        query = 'INSERT INTO public.contiguous_id_blocks VALUES ' + serialise_tuple(row)
        print(query)
        if args.apply_queries:
            prod_cursor.execute(query)

    # Get the last id used
    query = 'select max(id) from public.contiguous_id_blocks'
    prod_cursor.execute(query)
    last_id_used, = prod_cursor.fetchall()[0]

    # Set it to the sequence
    query = f"SELECT setval('hibernate_sequence', {last_id_used}, true);"
    print(query)
    if args.apply_queries:
        prod_cursor.execute(query)

    if args.apply_queries:
        prod_cursor.commit()
    prod_conn.close()
    dev_conn.close()


if __name__ == '__main__':
    main()
