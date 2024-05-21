import os

import requests
from requests.auth import HTTPBasicAuth

base_url = 'http://localhost:8081/eva/webservices/seqcol/'
admin_username = ''
admin_password = ''

files_to_download = [
    'https://raw.githubusercontent.com/refgenie/seqcolapi/master/analysis/data/test_data/base.fa',
    'https://raw.githubusercontent.com/refgenie/seqcolapi/master/analysis/data/test_data/different_names.fa',
    'https://raw.githubusercontent.com/refgenie/seqcolapi/master/analysis/data/test_data/different_order.fa',
    'https://raw.githubusercontent.com/refgenie/seqcolapi/master/analysis/data/test_data/pair_swap.fa',
    'https://raw.githubusercontent.com/refgenie/seqcolapi/master/analysis/data/test_data/subset.fa',
    'https://raw.githubusercontent.com/refgenie/seqcolapi/master/analysis/data/test_data/swap_wo_coords.fa'
]


def download_file(file_to_download):
    response = requests.get(file_to_download)
    response.raise_for_status()
    return os.path.basename(file_to_download), response.text


def upload_file(name, content):
    url = base_url + f'admin/seqcols/fasta/{name}'
    response = requests.put(url, data=content, auth=HTTPBasicAuth(admin_username, admin_password))
    response.raise_for_status()
    return response.json()


for f in files_to_download:
    name, content = download_file(f)
    result = upload_file(name, content.encode('ascii'))
    if result:
        # There should only be one
        print(result['seqcols'][0]['digest'])




