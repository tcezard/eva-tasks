import argparse
import csv
import getpass
from os import path

from ebi_eva_common_pyutils.logger import logging_config
from pymongo import MongoClient

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


def read_chunks_info(uri, mongo_password):
    mongo_client = MongoClient(uri, password=mongo_password)
    config_db = mongo_client["config"]
    chunk_collection = config_db["chunks"]

    skip_config_db = {
        "$match":
            {
                "ns": {
                    "$nin": ["config.system.sessions"],
                }
            }
    }
    group_by_collection_and_shard = {
        "$group":
            {
                "_id": {"collection": "$ns", "shard": "$shard"},
                "chunks": {"$sum": 1},
            }
    }
    sort_by_db_collection_shard = {
        "$sort":
            {
                "_id.collection": 1,
                "_id.shard": 1
            }
    }

    pipeline = [
        skip_config_db,
        group_by_collection_and_shard,
        sort_by_db_collection_shard
    ]

    logger.info("reading chunk info from config db's chunk collection...")
    aggregated_chunk_info = chunk_collection.aggregate(pipeline, batchSize=0)

    result = {}
    for agg_chunk in aggregated_chunk_info:
        ns = agg_chunk["_id"]["collection"]
        curr_shard = agg_chunk["_id"]["shard"]
        no_of_chunks = agg_chunk["chunks"]

        if ns in result.keys():
            result[ns]["shards"][curr_shard]["no_of_chunks"] = no_of_chunks
        else:
            database = ns.split(".")[0]
            collection = ns.split(".")[1]

            logger.info(f"reading collection info for db({database}) and collection({collection})")
            collection_info = mongo_client[database].command("collstats", collection)
            collection_size = int(collection_info["size"])
            collection_size_gb = round(collection_size / 1024 / 1024 / 1024, 2)

            result[ns] = {
                "database": database,
                "collection": collection,
                "collection_size": collection_size,
                "collection_size_gb": collection_size_gb,
                "shards": {}
            }

            shard_dict = collection_info["shards"]
            for shard_key in sorted(shard_dict.keys()):
                shard_size = int(shard_dict[shard_key]["size"])
                shard_size_gb = round(shard_size / 1024 / 1024 / 1024, 2)
                result[ns]["shards"][shard_key] = {
                    "shard_size": shard_size,
                    "shard_size_gb": shard_size_gb
                }

            result[ns]["shards"][curr_shard]["no_of_chunks"] = no_of_chunks

    mongo_client.close()

    return result


def write_chunk_info_to_csv(result_data, report_dir):
    rows = [["Database", "Collection", "Collection_Size - Bytes", "Collection_Size - GB", "Shard", "Shard_Size - Bytes",
             "Shard_Size - GB", "No_Of_Chunks"]]

    for key, value in result_data.items():
        database = value["database"]
        collection = value["collection"]
        collection_size = value["collection_size"]
        collection_size_gb = value["collection_size_gb"]
        for shard_key in sorted(value["shards"].keys()):
            shard = value["shards"][shard_key]
            print(database, collection, shard)
            rows.append([database, collection, collection_size, collection_size_gb, shard_key,
                         shard["shard_size"], shard["shard_size_gb"], shard["no_of_chunks"]])

    with open(path.join(report_dir, 'Shard_Report.csv'), 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description='Generate Reports on how collections of a db is sharded in a given MongoDB source',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--mongo-source-uri",
                        help="Mongo Source URI (ex: mongodb://user:@mongos-source-host:27017/admin)", required=True)
    parser.add_argument("--report-dir",
                        help="Top-level directory where report will be saved (ex: /path/to/report_dir)",
                        required=True)
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    mongo_password = getpass.getpass(prompt='Please Enter Mongo Source Password: ')
    result = read_chunks_info(args.mongo_source_uri, mongo_password)
    write_chunk_info_to_csv(result, args.report_dir)


if __name__ == "__main__":
    main()
