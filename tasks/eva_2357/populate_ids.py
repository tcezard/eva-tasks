import hashlib
import argparse
import pymongo
import traceback
import logging
from ebi_eva_common_pyutils.config_utils import get_mongo_uri_for_eva_profile


def generate_update_statement(hash_to_variant_ids, hash_to_accession_info):
    variant_to_ids = {}
    update_statements = []
    for sve_hash, accessions in hash_to_accession_info.items():
        ss_accession = f"ss{accessions.get('ss')}"
        rs_accession = f"rs{accessions.get('rs')}" if accessions.get('rs') else accessions.get('rs')
        variant_id = hash_to_variant_ids.get(sve_hash)

        # Group all the ids (ss and rs) by variant id
        if variant_to_ids.get(variant_id):
            variant_to_ids[variant_id].add(ss_accession)
            if rs_accession:
                variant_to_ids[variant_id].add(rs_accession)
        else:
            ids = set()
            ids.add(ss_accession)
            if rs_accession:
                ids.add(rs_accession)
            variant_to_ids[variant_id] = ids

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
    cursor_accessioning = sve_collection.find(sve_filter, projection=sve_projection)
    hash_to_accession_info = {}
    for sve in cursor_accessioning:
        sve_hash = sve.get("_id")
        ss_accession = sve.get('accession')
        rs_accession = sve.get('rs')
        hash_to_accession_info[sve_hash] = {"ss": ss_accession}
        # Get RS ID is variant is clustered
        if rs_accession:
            hash_to_accession_info[sve_hash]['rs'] = rs_accession
    return hash_to_accession_info


def get_variants_from_variant_warehouse(variants_collection):
    projection = {"_id": 1, "files.sid": 1, "chr": 1, "start": 1, "ref": 1, "alt": 1, "type": 1}
    return variants_collection.find({}, projection=projection)


def load_synonyms_for_assembly(assembly_accession, assembly_report_file=None):
    """
    Reads an assembly report and loads dictionaries to map names to genbank.
    Returns 5 dictionaries: by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc
    Example usage:
    genbank = by_name['1']['genbank']
    """
    logging.info(f"Parsing assembly report {assembly_report_file}")
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
            is_chromosome = columns[1] == 'assembled-molecule' and columns[3] == 'Chromosome'
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

    logging.info('Loaded chromosome synonyms for assembly {}'.format(assembly_accession))
    return by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc


def get_genbank(synonym_dictionaries, contig):
    """Returns a tuple (genbank, was_already_genbank) or raises an exception if the contig was not found"""
    by_name, by_assigned_molecule, by_genbank, by_refseq, by_ucsc = synonym_dictionaries
    if contig in by_genbank:
        return contig, True
    if contig in by_name:
        return by_name[contig]['genbank'], False
    if contig in by_assigned_molecule:
        return by_assigned_molecule[contig]['genbank'], False
    if contig in by_ucsc:
        return by_ucsc[contig]['genbank'], False
    if contig in by_refseq and by_refseq[contig]['is_genbank_refseq_identical']:
        return by_refseq[contig]['genbank'], False
    raise Exception(f"Could not find synonym for contig {contig}")


def get_SHA1(variant_id):
    """Calculate the SHA1 digest from the seq, study, contig, start, ref, and alt attributes of the variant"""
    return hashlib.sha1(variant_id.encode()).hexdigest().upper()


def get_ids(assembly, contig_synonym_dictionaries, variant_query_result, hash_to_variant_id):
    """
    Get a dictionary with the submitted variant hash as key and the variant id (chr_start_ref_alt) from the variant
    warehouse as value
    """
    for file in variant_query_result['files']:
        study = file['sid']
        start = variant_query_result['start']
        ref = variant_query_result['ref']
        alt = variant_query_result['alt']
        id_variant_warehouse = variant_query_result['_id']
        try:
            genbank_chr, was_already_genbank = get_genbank(contig_synonym_dictionaries, variant_query_result['chr'])
        except Exception as e:
            logging.info(f"{e} in variant {variant_query_result['chr']}_{start}_{ref}_{alt}")
            continue
        sve_hash = get_SHA1(f"{assembly}_{study}_{genbank_chr}_{start}_{ref}_{alt}")
        hash_to_variant_id[sve_hash] = id_variant_warehouse
    return hash_to_variant_id


def get_db_name_and_assembly_accession(databases):
    db_assembly = {}
    if databases:
        with open(databases, 'r') as db_file:
            for line in db_file:
                db, assembly, asm_report = line.split(',')
                db_assembly[db] = {"assembly": assembly, "asm_report": asm_report}
        return db_assembly


def populate_ids(private_config_xml_file, databases, profile='production', mongo_accession_db='eva_accession_sharded'):
    db_assembly = get_db_name_and_assembly_accession(databases)
    for db_name, info in db_assembly.items():
        assembly = info['assembly']
        asm_report = info['asm_report']
        logging.info(f"Processing database {db_name} (assembly {assembly})")

        contig_synonym_dictionaries = load_synonyms_for_assembly(assembly, asm_report)

        with pymongo.MongoClient(get_mongo_uri_for_eva_profile(profile, private_config_xml_file)) as mongo_handle:
            variants_collection = mongo_handle[db_name]["variants_2_0"]
            logging.info(f"Querying variants from variant warehouse, database {db_name}")
            variants_cursor = get_variants_from_variant_warehouse(variants_collection)
            hash_to_variant_id = {}
            update_statements = []
            try:
                for variant_query_result in variants_cursor:
                    # One variant in the variant warehouse can include more than one study, which means there can be
                    # more than one submitted variant per variant warehouse variant
                    get_ids(assembly, contig_synonym_dictionaries, variant_query_result, hash_to_variant_id)
                sve_hashes = hash_to_variant_id.keys()
                logging.info(f"Querying accessioning database for assembly {assembly}")
                hash_to_accession_info = get_from_accessioning_db(mongo_handle, mongo_accession_db, sve_hashes)
                update_statements.extend(generate_update_statement(hash_to_variant_id, hash_to_accession_info))
            except Exception as e:
                print(traceback.format_exc())
                raise e
            finally:
                variants_cursor.close()
        logging.info(f"Updating database {db_name} variant warehouse with ss and rs ids")
        result_update = variants_collection.bulk_write(requests=update_statements, ordered=False)
        logging.info(f"{result_update.modified_count} variants modified in {db_name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get stats from variant warehouse', add_help=False)
    parser.add_argument("--private-config-xml-file", help="ex: /path/to/eva-maven-settings.xml", required=True)
    parser.add_argument("--dbs-to-populate-list",
                        help="Full path to file with list of pairs (database, assembly) to populate the ids field "
                             "with SS and RS ids (ex: /path/to/dbs/to/populate.txt)", required=True)
    args = parser.parse_args()
    populate_ids(args.private_config_xml_file, args.dbs_to_populate_list)
