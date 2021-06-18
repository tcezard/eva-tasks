import hashlib
import argparse
from collections import defaultdict

import pymongo
import traceback
from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
from pymongo import WriteConcern

logging_config.add_stdout_handler()
logger = logging_config.get_logger(__name__)


def generate_update_statement(hash_to_variant_ids, hash_to_accession_info):
    variant_to_ids = defaultdict(set)
    update_statements = []

    # Group ss and rs accessions by variant id (chr_start_ref_alt)
    for sve_hash, accessions in hash_to_accession_info.items():
        variant_id = hash_to_variant_ids.get(sve_hash)
        variant_to_ids[variant_id].update(accessions)

    # Generate one update statement per variant_id
    for variant_id, ids in variant_to_ids.items():
        query = {"_id": variant_id}
        update = {"$addToSet": {"ids": {"$each": list(ids)}}}
        update_statements.append(pymongo.UpdateOne(query, update))
    return update_statements


def get_from_accessioning_db(mongo_handle, mongo_accession_db, sve_hashes):
    # Get SS ID from accessioning DB
    sve_collection = mongo_handle[mongo_accession_db]['submittedVariantEntity']
    sve_filter = {"_id": {"$in": list(sve_hashes)}}
    sve_projection = {"_id": 1, "accession": 1, "rs": 1}
    cursor_accessioning = sve_collection.find(sve_filter, projection=sve_projection, no_cursor_timeout=True)
    hash_to_accession_info = defaultdict(list)
    for sve in cursor_accessioning:
        sve_hash = sve.get("_id")
        ss_accession = sve.get('accession')
        rs_accession = sve.get('rs')
        hash_to_accession_info[sve_hash].append(f"ss{ss_accession}")
        # Get RS ID is variant is clustered
        if rs_accession:
            hash_to_accession_info[sve_hash].append(f"rs{rs_accession}")
    return hash_to_accession_info


def get_variants_from_variant_warehouse(variants_collection, batch_size):
    projection = {"_id": 1, "files.sid": 1, "chr": 1, "start": 1, "ref": 1, "alt": 1, "type": 1}
    return variants_collection.find({}, projection=projection, batch_size=batch_size, no_cursor_timeout=True)


def load_synonyms_for_assembly(assembly_accession, assembly_report_file=None):
    """
    Reads an assembly report and loads dictionaries to map names to genbank.
    Returns 5 dictionaries: by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc
    Example usage:
    genbank = by_name['1']['genbank']
    """
    logger.info(f"Parsing assembly report {assembly_report_file}")
    by_name = dict()
    by_assigned_molecule = dict()
    by_genbank = dict()
    by_refseq = dict()
    by_ucsc = dict()
    with open(assembly_report_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            is_chromosome = columns[1] == 'assembled-molecule' and columns[3].lower() == 'chromosome'
            synonyms = {'name': columns[0],
                        'assigned_molecule': columns[2] if is_chromosome else None,
                        'is_chromosome': is_chromosome,
                        'genbank': columns[4],
                        'refseq': columns[6],
                        'is_genbank_refseq_identical': (columns[5] == '='),
                        'ucsc': columns[9] if columns[9] != 'na' else None}

            by_name[columns[0]] = synonyms
            by_genbank[columns[4]] = synonyms
            by_refseq[columns[6]] = synonyms

            if synonyms['assigned_molecule'] is not None:
                by_assigned_molecule[columns[2]] = synonyms

            if synonyms['ucsc'] is not None:
                by_ucsc[columns[9]] = synonyms

    logger.info('Loaded chromosome synonyms for assembly {}'.format(assembly_accession))
    return by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc


def get_genbank(synonym_dictionaries, contig):
    """Returns the genbank accession or raises an exception if the contig was not found"""
    by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc = synonym_dictionaries
    if contig in by_genbank:
        return contig
    if contig in by_name:
        return by_name[contig]['genbank']
    if contig in by_assigned_molecule:
        return by_assigned_molecule[contig]['genbank']
    if contig in by_ucsc:
        return by_ucsc[contig]['genbank']
    if contig in by_refseq and by_refseq[contig]['is_genbank_refseq_identical']:
        return by_refseq[contig]['genbank']
    raise KeyError(f"Could not find synonym for contig {contig}")


def get_SHA1(variant_id):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    return hashlib.sha1(variant_id.encode()).hexdigest().upper()


def get_hash_to_variant_id(assembly, contig_synonym_dictionaries, variant_query_result):
    """
    Get a dictionary with the submitted variant hash as key and the variant id (chr_start_ref_alt) from the variant
    warehouse as value
    """
    hash_to_variant_id = {}
    contigs_no_genbank = []
    # One variant in the variant warehouse can include more than one study, which means there can be more than one
    # submitted variant per variant warehouse variant
    for file in variant_query_result['files']:
        study = file['sid']
        start = variant_query_result['start']
        ref = variant_query_result['ref']
        alt = variant_query_result['alt']
        id_variant_warehouse = variant_query_result['_id']
        chr = variant_query_result['chr']
        try:
            genbank_chr = get_genbank(contig_synonym_dictionaries, chr)
        except KeyError as e:
            contigs_no_genbank.append(chr)
            logger.error(f"{e} in variant {variant_query_result['chr']}_{start}_{ref}_{alt}")
            raise ValueError(f"Contig {chr} don't have a genbank equivalent, check assembly report")
        sve_hash = get_SHA1(f"{assembly}_{study}_{genbank_chr}_{start}_{ref}_{alt}")
        hash_to_variant_id[sve_hash] = id_variant_warehouse
    return hash_to_variant_id, contigs_no_genbank


def get_db_name_and_assembly_accession(databases):
    db_assembly = {}
    if databases:
        with open(databases, 'r') as db_file:
            for line in db_file:
                db, assembly, asm_report = line.split(',')
                db_assembly[db] = {"assembly": assembly, "asm_report": asm_report}
        return db_assembly


def update_variant_warehouse(mongo_handle, mongo_accession_db, variants_collection, hash_to_variant_id):
    update_statements = []
    sve_hashes = hash_to_variant_id.keys()
    hash_to_accession_info = get_from_accessioning_db(mongo_handle, mongo_accession_db, sve_hashes)
    update_statements.extend(generate_update_statement(hash_to_variant_id, hash_to_accession_info))
    if update_statements:
        result_update = variants_collection.with_options(
            write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
            .bulk_write(requests=update_statements, ordered=False)
        variants_modified_in_batch = result_update.modified_count if result_update else 0
        return variants_modified_in_batch


def populate_ids(private_config_xml_file, databases, profile='production', mongo_accession_db='eva_accession_sharded'):
    db_assembly = get_db_name_and_assembly_accession(databases)
    for db_name, info in db_assembly.items():
        assembly = info['assembly']
        asm_report = info['asm_report']
        logger.info(f"Processing database {db_name} (assembly {assembly})")

        contig_synonym_dictionaries = load_synonyms_for_assembly(assembly, asm_report)

        with pymongo.MongoClient(get_mongo_uri_for_eva_profile(profile, private_config_xml_file)) as mongo_handle:
            variants_collection = mongo_handle[db_name]["variants_2_0"]
            logger.info(f"Querying variants from variant warehouse, database {db_name}")
            batch_size = 1
            variants_cursor = get_variants_from_variant_warehouse(variants_collection, batch_size)
            hash_to_variant_ids = {}
            update_statements = []
            try:
                count_variants = 0
                batch_number = 0
                for variant_query_result in variants_cursor:
                    hash_to_variant_id, contigs_no_genbank = get_hash_to_variant_id(assembly,
                                                                                    contig_synonym_dictionaries,
                                                                                    variant_query_result)
                    hash_to_variant_ids.update(hash_to_variant_id)
                    count_variants += 1
                    if count_variants == batch_size:
                        batch_number += 1
                        # Generate update statements
                        logger.info(f"Generating update statements: database {db_name} (batch {batch_number})")
                        sve_hashes = hash_to_variant_ids.keys()
                        hash_to_accession_info = get_from_accessioning_db(mongo_handle, mongo_accession_db, sve_hashes)
                        update_statements.extend(generate_update_statement(hash_to_variant_ids, hash_to_accession_info))

                        hash_to_variant_ids.clear()
                        count_variants = 0
                if count_variants > 0:
                    logger.info(f"Generating update statements: database {db_name} (batch {batch_number+1})")
                    sve_hashes = hash_to_variant_ids.keys()
                    hash_to_accession_info = get_from_accessioning_db(mongo_handle, mongo_accession_db, sve_hashes)
                    update_statements.extend(generate_update_statement(hash_to_variant_ids, hash_to_accession_info))
            except ValueError as e:
                print(traceback.format_exc())
                raise e
            finally:
                variants_cursor.close()

            if len(contigs_no_genbank) > 0:
                raise ValueError(f"Contigs {contigs_no_genbank} don't have a genbank equivalent, check assembly report")

            modified_count = 0
            if update_statements:
                result_update = variants_collection.with_options(
                    write_concern=WriteConcern(w="majority", wtimeout=1200000)) \
                    .bulk_write(requests=update_statements, ordered=False)
                modified_count = result_update.modified_count if result_update else 0
                logger.info(f"{modified_count} variants modified in {db_name}")

        return modified_count


def check_all_contigs(private_config_xml_file, databases, profile='production'):
    db_assembly = get_db_name_and_assembly_accession(databases)
    for db_name, info in db_assembly.items():
        assembly = info['assembly']
        asm_report = info['asm_report']
        logger.info(f"Check database {db_name} (assembly {assembly}) with report {asm_report}")

        contig_synonym_dictionaries = load_synonyms_for_assembly(assembly, asm_report)

        with pymongo.MongoClient(get_mongo_uri_for_eva_profile(profile, private_config_xml_file)) as mongo_handle:
            variants_collection = mongo_handle[db_name]["variants_2_0"]
            cursor = variants_collection.aggregate([{'$group': {'_id': '$chr', 'count': {'$sum': 1}}}])

        translatable_contigs = 0
        translatable_variants = 0
        notranslation_variants = 0
        notranslation_contigs = set()
        for contig in cursor:
            try:
                _ = get_genbank(contig_synonym_dictionaries, contig['_id'])
                translatable_contigs += 1
                translatable_variants += contig['count']
            except KeyError:
                notranslation_variants += contig['count']
                notranslation_contigs.add(contig['_id'])

        if len(notranslation_contigs) > 0:
            raise ValueError(f'Aborting Update (no changes were done). '
                             f'With the provided assembly_report, the next {len(notranslation_contigs)} '
                             f'contigs (present in {notranslation_variants} variants) can not be '
                             f'replaced: {notranslation_contigs}')
        else:
            logger.info(f'Check ok. {translatable_contigs} contigs of {translatable_variants} variants will be  '
                        f'translated to genbank')
        return translatable_variants


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        add_help=False,
        description='Update variant warehouse with SS and RS id from the accessioning warehouse'
    )
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--dbs-to-populate-list",
                        help="Full path to file with list of pairs (database, assembly) to populate the ids field "
                             "with SS and RS ids (ex: /path/to/dbs/to/populate.txt)", required=True)
    parser.add_argument('--only_check', help='Check the contigs have a genbank equivalent in the assembly report',
                        default=False, action='store_true')
    parser.add_argument('--fail-on-first-error', help='Stop execution if one contig does not have a genbank equivalent',
                        default=False, action='store_true')
    args = parser.parse_args()

    check_all_contigs(args.private_config_xml_file, args.dbs_to_populate_list)
    if not args.only_check:
        populate_ids(args.private_config_xml_file, args.dbs_to_populate_list)
