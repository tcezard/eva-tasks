import requests

vcf_validator_releases = 'https://api.github.com/repos/EBIvariation/vcf-validator/releases'
response = requests.get(vcf_validator_releases)
response.raise_for_status()
releases = response.json()
for release in releases:
    release_id = release.get('id')
    response = requests.get(vcf_validator_releases + '/' + str(release_id))
    response.raise_for_status()
    release_detail = response.json()
    out = []
    for asset in release_detail.get('assets'):
        print('\t'.join([
            release_detail.get('tag_name'),
            asset.get('name'),
            str(asset.get('download_count'))
        ]))

