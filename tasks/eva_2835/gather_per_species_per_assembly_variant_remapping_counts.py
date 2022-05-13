#Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Importing relevant packages
import os
import yaml
import argparse
import xlsxwriter
from string import ascii_uppercase


def generate_output(all_columns, all_columns_for_specific_taxid_assembly, all_values_for_specific_taxid_assembly,
                    all_tax_id, all_assembly, output_path):
    """
    This function is used to generate the results in a spreadsheet containing the statistics of remapping
    and various reasons of remapping failures for different steps(i.e. flanking region length) for all the
    taxonomies and their corresponding assemblies

    Input: It accepts the remapping full path and the output path along with the Excel headers and its
    corresponding values

    Output: It generates an Excel spreadsheet with the statistics

    """

    # Initializing the Excel workbook
    workbook = xlsxwriter.Workbook(output_path + '/Gather_Stats.xlsx')
    worksheet = workbook.add_worksheet('All counts.xlsx')

    # Creating the header format to be used during writing of the results
    header_format = workbook.add_format({
        'bold': 1,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter',
        'fg_color': '#D7E4BC',
        'text_wrap': 1})

    # Considering maximum reasons of failure to be restricted to 26*3
    letters = list(ascii_uppercase[:3])
    all_letters = list(ascii_uppercase)
    excel_cell_headers = list(ascii_uppercase)

    # Storing the Excel cell header subscripts in a list
    for char in letters:
        for char2 in all_letters:
            excel_cell_headers.append(char + char2)

    # Displaying the first two headers i.e. the taxonomy id and the assembly accession
    worksheet.write('A1', 'Taxonomy', header_format)
    worksheet.write('B1', 'Assembly Accession', header_format)

    # Displaying the other column headers
    for column in range(len(all_columns)):
        cell_name = excel_cell_headers[column + 2] + "1"
        worksheet.write(cell_name, all_columns[column], header_format)

    # Defining the cell format for the taxonomy ids and the assembly accession columns
    filename_fmt = workbook.add_format(
        {'align': 'left', 'valign': 'vcenter', 'border': 1, 'font_color': 'red', 'text_wrap': 1})

    # Defining the cell format for the column values
    stats_fmt = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1})

    # Setting the starting row for entering the records
    row = 1

    # Populating the Excel with the values per taxonomy id and per assembly
    for r in range(len(all_tax_id)):
        worksheet.write(row, 0, all_tax_id[row - 1], filename_fmt)
        worksheet.write(row, 1, all_assembly[row - 1], filename_fmt)

        # Setting the index list which has values for a particular iteration
        idx = []

        for column in all_columns_for_specific_taxid_assembly[row - 1]:
            position = all_columns.index(column)
            idx.append(position)

        # Populating the values for the other stats columns for a particular iteration
        for col in range(len(all_columns)):
            if col in idx:
                key = all_columns[col]
                position = all_columns_for_specific_taxid_assembly[row - 1].index(key)
                stats = all_values_for_specific_taxid_assembly[row - 1][position]
                worksheet.write(row, col + 2, stats, stats_fmt)

            else:
                worksheet.write(row, col + 2, 0, stats_fmt)

        # Incrementing the row to write the next record
        row = row + 1

    # Closing the workbook
    workbook.close()


def extract_taxid_assembly(remapping_root_path):
    """
    This function is used to extract the full list of taxonomy ids and assembly accessions and pass it to
    gather_counts_per_tax_per_assembly() function for a particular taxonomy id and assembly for gathering statistics and subsequently consolidates statistics for a holistic analysis

    Input: It accepts the remapping full path for subsequent analysis

    Output: It returns the list of taxonomies and their corresponding assemblies in a dictionary format
    """

    # Collecting the tax_ids from the input path
    taxids = [name for name in os.listdir(remapping_root_path) if
              os.path.isdir(os.path.join(remapping_root_path, name))]

    # Defining a dictionary to store the taxonomy and its corresponding assemblies
    tax_assembly = {}

    for tax_id in taxids:
        assembly_accession = [name for name in os.listdir(os.path.join(remapping_root_path, tax_id)) if
                              os.path.isdir(os.path.join(os.path.join(remapping_root_path, tax_id), name))]
        tax_assembly[tax_id] = assembly_accession

    # Returning the dictionary containing the taxonomy and its corresponding assembly
    return tax_assembly


def gather_counts_per_tax_per_assembly(path, taxid, assembly_accession):
    """
    This function is used to store the counts of the remapped variants along wih the reason for failures
    at different rounds (i.e. for different lengths of the flanking region) for a particular taxonomy and
    assembly

    Input: The taxonomy id and the assembly accession from generate_output() along with the input files full
           path

    Output: All the statistics i.e. columns and values for a particular taxonomy and assembly

    """

    # Setting the filename for the eva counts
    filename_eva = str(assembly_accession) + "_eva_remapped_counts.yml"
    filename_eva = os.path.join(path, str(taxid), str(assembly_accession), "eva", filename_eva)

    # Setting the filename for the dbsnp counts
    filename_dbsnp = str(assembly_accession) + "_dbsnp_remapped_counts.yml"
    filename_dbsnp = os.path.join(path, str(taxid), str(assembly_accession), "dbsnp", filename_dbsnp)

    # Calling functions for gathering stats per file from eva and dbsnp data and storing the data in a dictionary
    # as a key-value pair
    columns_values_eva = gather_counts_per_file(filename_eva)
    columns_values_dbsnp = gather_counts_per_file(filename_dbsnp)

    # Adding the values of dbsnp and eva with the common columns i.e. fields and also concatenating them
    columns_values = {}
    for column in set(columns_values_eva) | set(columns_values_dbsnp):
        columns_values[column] = columns_values_eva.get(column, 0) + columns_values_dbsnp.get(column, 0)

    # Sorting the dictionary using the keys in order to maintain uniformity in the order of appearance of columns
    # in the final output spreadsheet
    columns_values = {column: value for column, value in sorted(columns_values.items())}

    # Splitting the dictionary into two lists - one for keys and one for values
    columns, values = zip(*columns_values.items())
    columns = list(columns)
    values = list(values)

    # Returning the fields in a yaml file along with its value for a particular tax id and accession
    return columns, values, taxid, assembly_accession


def gather_counts_per_file(filename):
    """
    This function is used to load the yml files for both eva and dbsnp and gather the data from the key-value pairs
    in the yml files

    Input: The absolute filepaths and filenames of the eva and the dbsnp files

    Output: A dictionary containing the key-value pairs of the yml files for a particular taxonomy and assembly

    """

    # Storing the contents for a particular yaml file in a linear dictionary format for eva

    # Defining a list to store all the fields in the yml files for a particular taxonomy and assembly
    columns = []

    # Defining a list to store all the values corresponding to the key fields in the yml files for a particular
    # taxonomy and assembly
    values = []

    with open(filename, 'r') as file:

        # Loading the data from the yaml file
        data = yaml.safe_load(file)

        # Iterating over the yml data for a particular file
        for column, value in data.items():

            # Renaming the column in order to sort the columns in proper order starting with Flank_050 data
            # followed by Flank_2000 data and finally by Flank_50000 data
            if column == "Flank_50":
                column = "Flank_050"

            # Checking if the value is dictionary
            if isinstance(value, dict):

                # Concatenating the outer field name i.e.(Flank_050, Flank_2000, Flank_50000) with the inner field
                # name i.e. the reasons of failures, remapping rate or flank_total
                for column2, value2 in value.items():
                    columns.append(column + "_" + column2)
                    values.append(value2)
            else:

                # Capitalizing the Taxonomy and the Assembly Accession field names in order to preserve the sorting
                # order of columns
                columns.append(column.capitalize())
                values.append(value)

    # Creating dictionary from two lists i.e. the columns and the values
    zip_iterator = zip(columns, values)
    columns_values = dict(zip_iterator)

    # Returning the fields and values for a particular yml file in a dictionary format
    return columns_values


def main():
    """
    This is the main function which accepts two arguments from the user, namely, remapping_root_path and
    the output path respectively.

    """

    parser = argparse.ArgumentParser(
        description='Collecting statistics per taxonomy per assembly for variant remapping')

    # Allowing the user to enter the full path for the taxonomy folders which contain
    # the corresponding old assemblies which in turn contain the yaml files containing
    # the remapped statistics
    parser.add_argument("--remapping_root_path", type=str,
                        help="Path where the remapping directories are present", required=True)

    # Taking the output path from the user where the visualization files and output csv
    # files are present
    parser.add_argument("--output_path", type=str,
                        help="Path to the output .", required=True)

    args = parser.parse_args()

    # Initializing the final output columns for the output spreadsheet. These comprise the total number of
    # variants,remapped variants, total number of variants for each flanking regions and the various reasons of
    # failures for each flanking regions
    all_columns = []

    # Initializing the list to store all counts of remapped variants, total number of variants for each
    # flanking regions and the various reasons of failures for each flanking regions and their corresponding
    # counts for all taxonomies and assemblies
    all_columns_for_specific_taxid_assembly = []
    all_values_for_specific_taxid_assembly = []

    # Initializing the list to store the values of all taxonomy ids and assembly accessions
    all_tax_id = []
    all_assembly = []

    # Collecting the taxonomy-assembly information
    tax_assembly = extract_taxid_assembly(args.remapping_root_path)

    # Collect statistics for each taxonomy and each assembly
    for taxid, assembly in tax_assembly.items():
        for val in range(len(assembly)):
            columns, values, taxid, assembly_accession = gather_counts_per_tax_per_assembly(args.remapping_root_path,
                                                                                            taxid, assembly[val])

            # Adding the fields and values to the master list for each taxonomy and assembly
            all_columns_for_specific_taxid_assembly.append(columns)
            all_values_for_specific_taxid_assembly.append(values)

            # Adding the taxonomy and assembly accession for a particular file to the master list
            all_tax_id.append(taxid)
            all_assembly.append(assembly_accession)

            # Updating the final output column in a sorted order every time from an individual yml file
            all_columns.extend(columns)
            all_columns = list(dict.fromkeys(all_columns))
            all_columns.sort()

    # Generating the final output spreadsheet using all the statistics gathered for all taxonomies and assemblies
    generate_output(all_columns, all_columns_for_specific_taxid_assembly, all_values_for_specific_taxid_assembly,
                    all_tax_id, all_assembly, args.output_path)


if __name__ == "__main__":
    main()
