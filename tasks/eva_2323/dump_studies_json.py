import json
from argparse import ArgumentParser
from cached_property import cached_property

import psycopg2
from ebi_eva_common_pyutils.config_utils import get_properties_from_xml_file
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query


class StudyDumper:

    def __init__(self, properties_file, stage='development'):
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

    def dump(self):

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

    def get_fields(self, all_fields):
        return [
            {'name': f, 'value': all_fields[f]} for f in all_fields
        ]

    def get_entry(self, fields, cross_references, hierarchical_field):
        return {
            'fields': fields,
            'cross_references': cross_references,
            'hierarchical_field': hierarchical_field
        }

    def json_document(self, entries):
        return {
            'name': 'EVA',
            'entry_count': len(entries),
            'entries': entries
        }


def main():
    parser = ArgumentParser()
    parser.add_argument('--property_file')
    args = parser.parse_args()

    print(StudyDumper(args.property_file).dump())


if __name__ == "__main__":
    main()
