import gzip
from argparse import ArgumentParser

import requests


def get_rs_from_ss(ssid):
    response = requests.get(f'https://www.ebi.ac.uk/eva/webservices/identifiers/v1/submitted-variants/{ssid}')
    response.raise_for_status()
    json_data = response.json()
    assert len(json_data) == 1
    return json_data[0]['data'].get('clusteredVariantAccession') or json_data[0]['data'].get('backPropagatedVariantAccession')


parser = ArgumentParser()
parser.add_argument('--accessioning_report')
parser.add_argument('--annotated_vcf')
args = parser.parse_args()

input_vcf = args.accessioning_report
output_vcf = args.annotated_vcf
with open(output_vcf, 'w') as open_output:
    with gzip.open(input_vcf, 'rt') as open_file:
        for line in open_file:
            if line.startswith('#'):
                open_output.write(line)
            else:
                sp_line = line.split('\t')
                ssid = int(sp_line[2][2:])
                rsid = get_rs_from_ss(ssid)
                sp_line[2] = f'rs{rsid}'
                open_output.write('\t'.join(sp_line))
