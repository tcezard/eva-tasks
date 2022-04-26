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
import pandas as pd


# Initializing the final output column
all_columns=[]

# Initializing the list to store all the keys and values for all the files
all_keys=[]
all_values=[]

# Initializing the list to store the values of all taxonomy ids and assembly accessions
all_tax_id=[]
all_assembly=[]


def generate_output(remapping_root_path, output_path):
    """
    This function is used to store the full list of taxonomy ids and assembly accessions along with generation
    of the output excel containing all the stats gathered from the yml files

    Input: It accepts the remapping full path from the user along with the output path
           for storing the results for subsequent analysis

    Output: It calls gather_counts_per_tax_per_assembly() function
           for a particular taxonomy id and a particular assembly of that taxonomy and generate a spreadsheet
           containing the statistics of remapping and various reasons of remapping failures for different steps
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

    # Collect statistics for each taxonomy and each assembly
    for taxid, assembly in tax_assembly.items():
        for val in range(len(assembly)):
            gather_counts_per_tax_per_assembly(remapping_root_path, taxid, assembly[val])

    # Generate output file from the statistics gathered

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
    a_z = list(ascii_uppercase)

    # Storing the Excel cell header subscripts for Excel
    for char in letters:
        for char2 in all_letters:
            a_z.append(char + char2)

    # Displaying the first two headers
    worksheet.write('A1', 'Taxonomy', header_format)
    worksheet.write('B1', 'Assembly Accession', header_format)

    # Displaying the other headers
    for column in range(len(all_columns)):
        cell_name = a_z[column + 2] + "1"
        worksheet.write(cell_name, all_columns[column], header_format)

    # Defining the cell format for the values
    stats_fmt = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1})
    filename_fmt = workbook.add_format(
        {'align': 'left', 'valign': 'vcenter', 'border': 1, 'font_color': 'red', 'text_wrap': 1})

    # Setting the starting row for values
    row = 1

    # Populating the Excel with the values per taxonomy id and per assembly
    for r in range(len(all_tax_id)):
        worksheet.write(row, 0, all_tax_id[row - 1], filename_fmt)
        worksheet.write(row, 1, all_assembly[row - 1], filename_fmt)

        # Setting the index list which has values for a particular iteration
        idx = []

        for key in all_keys[row - 1]:
            position = all_columns.index(key)
            idx.append(position)

        # Populating the values for the other stats columns for a particular iteration
        for col in range(len(all_columns)):
            if col in idx:
                key = all_columns[col]
                position = all_keys[row - 1].index(key)
                stats = all_values[row - 1][position]
                worksheet.write(row, col + 2, stats, stats_fmt)

            else:
                worksheet.write(row, col + 2, 0, stats_fmt)

        # Incrementing the row to write the next record
        row = row + 1

    # Closing the workbook
    workbook.close()


def gather_counts_per_tax_per_assembly(path, taxid, assembly_accession):
    """
    This function is used to store the counts of the remapped variants along wih the reason for failures
    at different rounds (i.e. for different lengths of the flanking region) for a particular taxonomy and
    assembly

    Input: The taxonomy id and the assembly accession from generate_output() along with the input files full
           path

    Output: All the statistics in global lists

    """

    # Setting the filename for the eva counts
    filename_eva = str(assembly_accession) + "_eva_remapped_counts.yml"
    filename_eva = os.path.join(path, str(taxid), str(assembly_accession), "eva", filename_eva)

    # Setting the filename for the dbsnp counts
    filename_dbsnp = str(assembly_accession) + "_dbsnp_remapped_counts.yml"
    filename_dbsnp = os.path.join(path, str(taxid), str(assembly_accession), "dbsnp", filename_dbsnp)

    # Calling functions for gathering stats per file from eva and dbsnp data
    keys_values = gather_counts_per_file(filename_eva)
    keys_values_dbsnp = gather_counts_per_file(filename_dbsnp)

    # Adding the values of dbsnp and eva with the common keys
    for key in keys_values_dbsnp:
        if key in keys_values:
            keys_values_dbsnp[key] = keys_values_dbsnp[key] + keys_values[key]

    # Concatenating the two dictionaries
    keys_values_dbsnp = {**keys_values, **keys_values_dbsnp}

    # Sorting the dictionary using the keys
    keys_values_dbsnp = {key: value for key, value in sorted(keys_values_dbsnp.items())}

    # Splitting the dictionary into two lists - one for keys and one for values
    keys, values = zip(*keys_values_dbsnp.items())
    keys = list(keys)
    values = list(values)

    # Adding the keys and values list to the master list for each taxonomy and assembly
    all_keys.append(keys)
    all_values.append(values)

    # Adding the taxonomy and assembly accession for a particular file
    all_tax_id.append(taxid)
    all_assembly.append(assembly_accession)

    # Updating the final output columns
    global all_columns
    all_columns.extend(keys)
    all_columns = list(dict.fromkeys(all_columns)
    all_columns.sort()


def gather_counts_per_file(filename):
    """
    This function is used to load the yml files for both eva and dbsnp and gather the data for the kay-value pairs
    in the yml files

    Input: The absolute filepaths and filenames of the eva and the dbsnp files

    Output: A dictionary containing the key-value pairs of the yml files for a particular taxonomy and assembly

    """

    # Storing the contents for a particular yaml file in a linear dictionary format for eva

    # Defining a list to store all the fields in the yml files for a particular taxonomy and assembly
    keys = []

    # Defining a list to store all the values corresponding to the key fields in the yml files for a particular
    # taxonomy and assembly
    values = []

    # Defining a dictionary to store the key-value pairs of the yml files for a particular taxonomy and assembly
    keys_values = {}

    with open(filename, 'r') as file:

        # Loading the data from the yaml file
        data = yaml.safe_load(file)

        for key, value in data.items():

            if key == "Flank_50":
                key = "Flank_050"

            # Checking if the value is dictionary
            if isinstance(value, dict):
                for key2, value2 in value.items():
                    keys.append(key + "_" + key2)
                    values.append(value2)
            else:
                keys.append(key.capitalize())
                values.append(value)

    # Creating dictionary from two lists
    zip_iterator = zip(keys, values)
    keys_values = dict(zip_iterator)

    return keys_values


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

    # Calling the primary function responsible for generating the statistics
    generate_output(args.remapping_root_path, args.output_path)

if __name__ == "__main__":
    main()

