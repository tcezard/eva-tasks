# Copyright 2022 EMBL - European Bioinformatics Institute
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

import os
import argparse

def collect_structural_variant_types(gvf_root_path):
    all_variant_list = set()
    file_folder_list = os.listdir(gvf_root_path)
    file_list = [gvf_file for gvf_file in file_folder_list if os.path.isfile(os.path.join(gvf_root_path, gvf_file))]
    for file in file_list:
        if file.endswith(".gvf"):
            variant_list = []
            with open(os.path.join(gvf_root_path, file)) as open_file:
                for file_line in open_file:
                    if file_line.startswith("#"):
                        continue
                    attributes_column = file_line.split("\t")[2]
                    variant_list.append(attributes_column)

            variant_list_unique = set(variant_list)
            all_variant_list.update(variant_list_unique)
    return all_variant_list


def main():
    parser = argparse.ArgumentParser(
        description='Collecting different variant types present in the gvf files')

    parser.add_argument("--gvf_root_path", type=str,
                        help="Path where the gvf files are present", required=True)

    parser.add_argument("--output_path", type=str,
                        help="Path to the output .", required=True)
    args = parser.parse_args()
    variant_types = collect_structural_variant_types(args.gvf_root_path)
    with open(os.path.join(args.output_path, "output_file.txt"), 'w') as output_file:
        print("\nThe different variant types are as follows:", file=output_file)
        print(*variant_types, sep="\n", file=output_file)


if __name__ == "__main__":
    main()
