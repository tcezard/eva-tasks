import bz2
import click
import json
import traceback

from ebi_eva_common_pyutils.logger import logging_config

logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)


def is_rs_id_mapped_to_assembly(rs_record, eva_production_human_dbsnp_assembly):
    try:
        for placement in rs_record["primary_snapshot_data"]["placements_with_allele"]:
            # Primary top-level placement (ptlp) is used to identify records that are mapped to this assembly
            # See https://github.com/EBIvariation/eva-accession/blob/d3c849f15b6b63ebd97d4ac6f2d857e3b65d908f/eva-accession-import-dbsnp2/src/main/java/uk/ac/ebi/eva/accession/dbsnp2/batch/processors/JsonNodeToClusteredVariantProcessor.java#L74-L74
            if placement["is_ptlp"]:
                for assembly_info in placement["placement_annot"]["seq_id_traits_by_assembly"]:
                    if assembly_info["assembly_accession"] == eva_production_human_dbsnp_assembly:
                        return True
    except KeyError as ex:
        logger.error(traceback.format_exc())
        return False
    return False


def check_subsnp_assumptions(release_json_file, eva_production_human_dbsnp_assembly):
    with bz2.open(release_json_file) as release_json_file_handle:
        line_index = 0
        total_num_subsnps = 0
        for json_line in release_json_file_handle:
            if line_index % 100000 == 0:
                logger.info("Processed {0} records...".format(line_index))
            line_index += 1
            rs_record = json.loads(json_line.decode("utf-8").strip())
            if not (is_rs_id_mapped_to_assembly(rs_record, eva_production_human_dbsnp_assembly)):
                logger.error("RS ID {0} is not mapped to assembly {1}".format(rs_record["refsnp_id"],
                                                                              eva_production_human_dbsnp_assembly))
                continue
            rs_id = int(rs_record["refsnp_id"])
            all_subsnps_for_rs = []

            try:
                for movement_record in rs_record["present_obs_movements"]:
                    current_subsnps = []
                    if "component_ids" in movement_record:
                        for component_id_record in movement_record["component_ids"]:
                            if "type" in component_id_record:
                                if component_id_record["type"] == "subsnp":
                                    subsnp = int(component_id_record["value"])
                                    current_subsnps.append(subsnp)
                                    all_subsnps_for_rs.append(subsnp)
                    if len(all_subsnps_for_rs) > 0:
                        if "allele_in_cur_release" in movement_record:
                            allele_info = movement_record["allele_in_cur_release"]
                            if "seq_id" not in allele_info or "position" not in allele_info \
                                or "deleted_sequence" not in allele_info or "inserted_sequence" not in allele_info:
                                logger.info("{0} SS IDs: {1} do not have mapping information"
                                            .format(len(current_subsnps), ",".join([str(snp) for snp in current_subsnps]
                                                                                   )))
                        else:
                            logger.info("{0} SS IDs: {1} do not have mapping information"
                                        .format(len(current_subsnps), ",".join([str(snp) for snp in current_subsnps])))

            except KeyError as ex:
                pass

            if len(all_subsnps_for_rs) == 0:
                logger.error("No SS IDs found for RS ID {0}".format(rs_id))

            total_num_subsnps += len(set(all_subsnps_for_rs))

        logger.info("Total number of SubSNPs: {0}".format(total_num_subsnps))


@click.option("--release-json-file", help="ex: /path/to/release/json_file.json", required=True)
@click.option("--eva-production-human-dbsnp-assembly",
              help="Most recent dbSNP human assembly in EVA production (ex: GCF_000001405.38)", required=True)
@click.command()
def main(release_json_file, eva_production_human_dbsnp_assembly):
    check_subsnp_assumptions(release_json_file, eva_production_human_dbsnp_assembly)


if __name__ == "__main__":
    main()
