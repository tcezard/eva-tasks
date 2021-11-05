import argparse

from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle

from tasks.eva_2641.constants import annotation_metadata_collection_name


def get_affected_dbs(settings_xml_file):
    affected_dbs = []
    with get_mongo_connection_handle('production', settings_xml_file) as mongo_conn:
        for db_name in mongo_conn.list_database_names():
            if not db_name.startswith('eva_'):
                continue
            coll = mongo_conn[db_name][annotation_metadata_collection_name]
            count = coll.count_documents({'ct': {'$exists': True}})
            if count > 0:
                affected_dbs.append((db_name, count))
    return affected_dbs


def main():
    parser = argparse.ArgumentParser(
        description='Get list of databases with annotations in the wrong collection', add_help=False)
    parser.add_argument('--settings-xml-file', required=True)

    args = parser.parse_args()
    dbs = get_affected_dbs(args.settings_xml_file)
    for db, count in dbs:
        print(db, ':', count)


if __name__ == '__main__':
    main()
