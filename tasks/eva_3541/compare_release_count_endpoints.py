import requests

v1_per_assembly = 'https://www.ebi.ac.uk/eva/webservices/release/v1/stats/per-assembly?releaseVersion={version}'
v1_per_species = 'https://www.ebi.ac.uk/eva/webservices/release/v1/stats/per-species?releaseVersion={version}'
v2_per_assembly = 'https://wwwdev.ebi.ac.uk/eva/webservices/release/v2/stats/per-assembly?releaseVersion={version}'
v2_per_species = 'https://wwwdev.ebi.ac.uk/eva/webservices/release/v2/stats/per-species?releaseVersion={version}'


def compare_assembly_release_version(version):
    response = requests.get(v1_per_assembly.format(version=version))
    response.raise_for_status()
    assembly_data_v1 = response.json()
    response = requests.get(v2_per_assembly.format(version=version))
    response.raise_for_status()
    assembly_data_v2 = response.json()
    assembly_dict_v1 = dict([(a.get('assemblyAccession'), a) for a in assembly_data_v1])
    assembly_dict_v2 = dict([(a.get('assemblyAccession'), a) for a in assembly_data_v2])
    assemblies = set(assembly_dict_v1) | set(assembly_dict_v2)
    for assembly_accession in assemblies:
        assembly_data_v1 = assembly_dict_v1.get(assembly_accession)
        assembly_data_v2 = assembly_dict_v2.get(assembly_accession)
        if not assembly_data_v1:
            print(f'For release {version}, Assembly {assembly_accession} is missing in version 1 of the endpoint ')
            continue
        if not assembly_data_v2:
            print(f'For release {version}, Assembly {assembly_accession} is missing in version 2 of the endpoint ')
            continue
        for metric in ['currentRs', 'mergedRs', 'deprecatedRs', 'mergedDeprecatedRs']:
            out = [
                version, assembly_accession, metric,
                assembly_data_v1.get(metric), assembly_data_v2.get(metric),
                assembly_data_v1.get(metric) - assembly_data_v2.get(metric)
            ]
            diff_metric = assembly_data_v1.get(metric) - assembly_data_v2.get(metric)
            if diff_metric:
                # print('\t'.join([str(s) for s in out]))
                yield '\t'.join([str(s) for s in out])


def compare_species_release_version(version):
    response = requests.get(v1_per_species.format(version=version))
    response.raise_for_status()
    species_data_v1 = response.json()
    response = requests.get(v2_per_species.format(version=version))
    response.raise_for_status()
    species_data_v2 = response.json()
    species_dict_v1 = dict([(a.get('taxonomyId'), a) for a in species_data_v1] )
    species_dict_v2 = dict([(a.get('taxonomyId'), a) for a in species_data_v2] )
    taxonomies = set(species_dict_v1) | set(species_dict_v2)
    for taxonomy in taxonomies:
        scientific_name = None
        species_data_v1 = species_dict_v1.get(taxonomy)
        species_data_v2 = species_dict_v2.get(taxonomy)
        if species_data_v2:
            scientific_name = species_data_v2.get('scientificName')
        elif species_data_v1:
            scientific_name = species_data_v1.get('scientificName')
        if not species_data_v1:
            print(f'For release {version}, species {taxonomy} - {scientific_name} is missing in version 1 of the endpoint ')
            continue
        if not species_data_v2:
            print(f'For release {version}, species {taxonomy} - {scientific_name} is missing in version 2 of the endpoint ')
            continue

        for metric in ['currentRs', 'mergedRs', 'deprecatedRs', 'mergedDeprecatedRs']:
            out = [
                version, taxonomy, scientific_name, metric,
                species_data_v1.get(metric), species_data_v2.get(metric),
                species_data_v1.get(metric) - species_data_v2.get(metric)
            ]
            diff_metric = species_data_v1.get(metric) - species_data_v2.get(metric)
            if diff_metric:
                # print('\t'.join([str(s) for s in out]))
                yield '\t'.join([str(s) for s in out])




def check_all_versions():
    output_file = 'different_assembly_metrics.tsv'
    with open(output_file, 'w') as open_output:
        open_output.write('\t'.join(['Version', 'Assembly', 'Metric', 'Count v1', 'Count v2', 'Difference']) + '\n')
        for version in range(1, 6):
            different_metrics = compare_assembly_release_version(version)
            for line in different_metrics:
                open_output.write(line + '\n')

    output_file = 'different_species_metrics.tsv'
    with open(output_file, 'w') as open_output:
        open_output.write('\t'.join(['Version', 'Taxonomy', 'Scientific name', 'Metric', 'Count v1', 'Count v2', 'Difference']) + '\n')
        for version in range(1, 6):
            different_metrics = compare_species_release_version(version)
            for line in different_metrics:
                open_output.write(line + '\n')


if __name__ == '__main__':
    check_all_versions()
