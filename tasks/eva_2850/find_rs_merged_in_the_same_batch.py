import argparse
import os

from ebi_eva_common_pyutils.logger import logging_config

logger = logging_config.get_logger(__name__)


def merged_rs_ids_present_in_same_batch(all_log_files):
    merged_rs_ids = []
    merged_rs_ids_no_further_hash_collision = []
    merged_rs_with_further_hash_collision = []

    for log_file_path in sorted(all_log_files, key=lambda x: os.stat(x).st_size):
        with open(log_file_path, 'r') as log_file:
            curr_rs_id_batch = []
            for line in log_file:
                # if "Started processing assembly :" in line:
                #     logger.info(line)
                if "RS_list ->" in line:
                    # logger.info(curr_rs_id_batch)
                    rs_id_list = line[line.find('['):].split(",")
                    rs_id_list = [s.replace("[", "") for s in rs_id_list]
                    rs_id_list = [s.replace("]", "") for s in rs_id_list]
                    rs_id_list = [s.strip() for s in rs_id_list]
                    # only one rs in batch
                    if len(rs_id_list) == 1:
                        continue
                    curr_rs_id_batch = rs_id_list
                if "creating merge event for" in line:
                    # logger.info(line)
                    merged_rs = line[line.index("creating merge event for"):].split(":")[1].split(" ")[1].strip()
                    merged_into = line[line.index("creating merge event for"):].split(":")[2].strip()

                    # check if the variants (which should have been deleted) are involved in another hash collision
                    if merged_rs in merged_rs_ids:
                        merged_rs_with_further_hash_collision.append(merged_rs)
                    if merged_into in merged_rs_ids:
                        merged_rs_with_further_hash_collision.append(merged_into)

                    if merged_rs in curr_rs_id_batch and merged_into in curr_rs_id_batch:
                        # logger.info(curr_rs_id_batch)
                        if curr_rs_id_batch.index(merged_rs) > curr_rs_id_batch.index(merged_into):
                            logger.error(f"Merged RS : {merged_rs} Merged Into: {merged_into}")
                            # logger.error(f"grep -n -E \"{merged_rs}|{merged_into}\" {os.path.basename(log_file_path)}")
                            merged_rs_ids.append(merged_rs)

                if "No hash collision for RS" in line:
                    rs_id = line.split(" ")[-1].strip()
                    if rs_id in merged_rs_ids:
                        merged_rs_ids_no_further_hash_collision.append(rs_id)

    logger.info(f"All rs ids merged in same batch: {sorted(merged_rs_ids)}")
    logger.info(f"RS ids not involved in further merge events: "
                f"{sorted(set(merged_rs_ids_no_further_hash_collision) - set(merged_rs_with_further_hash_collision))}")
    logger.info(f"RS ids involved in further merge events: "
                f"{sorted(set(merged_rs_with_further_hash_collision) - set(merged_rs_ids_no_further_hash_collision))}")
    logger.info(f"RS ids processed twice - "
                f"once not involved in any merge and again involved in another merge event"
                f"{sorted(set(merged_rs_ids_no_further_hash_collision).intersection(set(merged_rs_with_further_hash_collision)))}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find variants which are part of the same batch and involved in merge',
                                     add_help=False)
    parser.add_argument("--log-file-dir", help="File containing discordant rs ids", required=True)
    args = parser.parse_args()

    # there are 2 different log files
    all_log_files = [os.path.join(args.log_file_dir, filename) for filename in os.listdir(args.log_file_dir)]

    merged_rs_ids_present_in_same_batch(all_log_files)
