import argparse
from collections import defaultdict

from ebi_eva_common_pyutils.mongodb import MongoDatabase


def categorise_many_ss_for_missing_rs(mongo_db, rsids, source, assembly_accession):
    if source == 'EVA':
        ss_source = 'submittedVariantEntity'
        op_source = 'clusteredVariantOperationEntity'
    elif source == 'DBSNP':
        ss_source = 'dbsnpSubmittedVariantEntity'
        op_source = 'dbsnpClusteredVariantOperationEntity'
    sve_collection = mongo_db.mongo_handle[mongo_db.db_name][ss_source]
    cvoe_collection = mongo_db.mongo_handle[mongo_db.db_name][op_source]
    cursor = sve_collection.find({'seq': assembly_accession, 'rs': {'$in': rsids}})
    rsid_to_type_set = defaultdict(set)
    for variant_rec in cursor:
        if 'remappedFrom' in variant_rec:
            rsid_to_type_set[variant_rec['rs']].add('Remapped cluster')
        else:
            rsid_to_type_set[variant_rec['rs']].add('Novo cluster')
    cursor.close()
    cursor = cvoe_collection.find({'inactiveObjects.asm': assembly_accession, 'accession': {'$in': rsids}})
    for variant_op in cursor:
        rsid_to_type_set[variant_op['accession']].add(variant_op.get('eventType'))
    cursor.close()
    cursor = cvoe_collection.find({'inactiveObjects.asm': assembly_accession, 'mergeInto': {'$in': rsids}})
    for variant_op in cursor:
        rsid_to_type_set[variant_op['mergeInto']].add('Target Of ' + variant_op.get('eventType'))
    cursor.close()
    cursor = cvoe_collection.find({'inactiveObjects.asm': assembly_accession, 'splitInto': {'$in': rsids}})
    for variant_op in cursor:
        rsid_to_type_set[variant_op['splitInto']].add('Target Of ' + variant_op.get('eventType'))
    cursor.close()
    rsid_to_types = {}
    for rsid in rsid_to_type_set:
        rsid_to_types[rsid] = sorted(rsid_to_type_set[rsid])
    return rsid_to_types


def categorise_all_rs(mongo_db, missing_rs_file, output_file):
    all_rs_accessions = defaultdict(set)
    with open(missing_rs_file) as open_file:
        for line in open_file:
            rs_accession, source, assembly = line.strip().split()
            all_rs_accessions[(source, assembly)].add(int(rs_accession))
    nb_processed = 0
    with open(output_file, 'w') as open_file:
        for source, assembly in all_rs_accessions:
            rsids = list(all_rs_accessions[(source, assembly)])
            rsid_to_types = categorise_many_ss_for_missing_rs(mongo_db, rsids, source, assembly)
            for rsid in rsid_to_types:
                print(f"{rsid}\t{source}\t{assembly}\t{','.join(rsid_to_types[rsid])}", file=open_file)
            nb_processed += len(rsids)
            print(f"{nb_processed}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check if the RS entities referenced in SS entities all exist in '
                                                 'the database')
    parser.add_argument("--missing_rs_file",
                        help="Input file containing the missing RSids and if it has ss in EVA or dbSNP", required=True)
    parser.add_argument("--output_file",
                        help="Output file containing the missing RSids annotated with the type", required=True)
    parser.add_argument("--mongo-db-uri",
                        help="Mongo Database URI (ex: mongodb://user:@mongos-db-host:27017/admin)", required=True)
    parser.add_argument("--mongo-db-secrets-file",
                        help="Full path to the Mongo Database secrets file (ex: /path/to/mongo/db/secret)",
                        required=True)

    args = parser.parse_args()
    mongo_db = MongoDatabase(uri=args.mongo_db_uri, secrets_file=args.mongo_db_secrets_file,
                             db_name="eva_accession_sharded")
    categorise_all_rs(mongo_db, args.missing_rs_file, args.output_file)
