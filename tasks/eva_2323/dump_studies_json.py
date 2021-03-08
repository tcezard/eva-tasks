import json
from argparse import ArgumentParser

import requests
from cached_property import cached_property

import psycopg2
from ebi_eva_common_pyutils.config_utils import get_properties_from_xml_file
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query


class StudyDumper:

    def __init__(self, properties_file=None, stage='production'):
        self.properties_file = properties_file
        self.stage = stage

    @cached_property
    def connect(self):
        properties = get_properties_from_xml_file(self.stage, self.properties_file)
        pg_url = properties['eva.evapro.jdbc.url'].split('jdbc:')[-1]
        pg_user = properties['eva.evapro.user']
        pg_pass = properties['eva.evapro.password']
        conn = psycopg2.connect(pg_url, user=pg_user, password=pg_pass)
        return conn

    def dump_from_database(self):

        query = (
            'SELECT p.project_accession, p.title, p.description, a.vcf_reference_accession, sa.taxonomy_id FROM project p '
            'JOIN project_analysis pa ON p.project_accession=pa.project_accession '
            'JOIN analysis a ON pa.analysis_accession=a.analysis_accession '
            'JOIN assembly_set sa ON sa.assembly_set_id=a.assembly_set_id'
        )
        rows = get_all_results_for_query(self.connect, query)
        entries = []
        for project_accession, title, description, vcf_reference_accession, taxonomy_id in rows:
            fields = self.get_fields({
                'id': project_accession,
                'title': title,
                'description': description,
                'genome_accession': vcf_reference_accession,
                'taxonomy_id': taxonomy_id
            })
            entries.append({'fields': fields})
        return json.dumps(self.json_document(entries), indent=4)

    def dump_from_api(self):
        url = 'https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/studies/all'
        response = requests.get(url)
        data = response.json()
        entries = []

        for study in data['response'][0]['result']:
            fields = self.get_fields({
                'id': study['id'],
                'title': study['name'],
                'description': study['description'],
                'genome_accession': study['assemblyAccession'],
                'taxonomy_id': study['taxonomyId'][0],
                'scientific_name': study['speciesScientificName'],
                'number_samples': study['numSamples']
            })
            cross_reference = []
            for publication in study['publications']:
                if publication and publication != '-':
                    cross_reference.append({'dbname': 'PUBMED', 'dbkey': str(publication)})
            cross_reference.append({'dbname': 'ENA', 'dbkey': study['id']})
            entries.append(self.get_entry(fields, cross_reference))
        return json.dumps(self.json_document(entries), indent=4)

    def get_fields(self, all_fields):
        return [
            {'name': f, 'value': all_fields[f]} for f in all_fields
        ]

    def get_entry(self, fields, cross_references=None, hierarchical_field=None):
        entry = {}
        if fields:
            entry['fields'] = fields
        if cross_references:
            entry['cross_references'] = cross_references
        if hierarchical_field:
            entry['hierarchical_field'] = hierarchical_field
        return entry

    def json_document(self, entries):
        return {
            'name': 'EVA studies',
            'entry_count': len(entries),
            'entries': entries
        }


def main():
    parser = ArgumentParser()
    parser.add_argument('--property_file', help='The properties file containing the credential to connect to EVAPRO. '
                                                'Without this option the REST API will be used')
    parser.add_argument('--stage',  choices=['development', 'production'], default='production',
                        help='The development of production server to be used for the database')
    args = parser.parse_args()
    if args.property_file:
        print(StudyDumper(args.property_file, args.stage).dump_from_database())
    else:
        print(StudyDumper().dump_from_api())


if __name__ == "__main__":
    main()
