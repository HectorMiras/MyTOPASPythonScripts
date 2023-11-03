from topas_csv_files_manager import merge_csv
from files_and_directory_manager import get_outputfile_paths
from phsp_manager import merge_phsp_files, merge_header_files
import json
import os
import sys

sim_path = sys.argv[1]
sim_config_file = sys.argv[2]
#sim_path = "/home/hector/mytopassimulations/MGHsimulations/tests/nodes_output_I125"
#sim_config_file = "/home/hector/mytopassimulations/MGHsimulations/tests/nodes_output_I125/simconfig.json"
output_dir = os.path.join(sim_path, "results")

#file_patterns = ["DoseToCell*", "DoseToNucleus*electrons.csv", "DoseToNucleus*gammas.csv", "nucleus_PHSP*"]

# If simconfig.json file exists in the directory, get the OUTPUT_FILE_PATTERNS from the file
if os.path.isfile(sim_config_file):
    # Open the file and read the JSON content
    with open(sim_config_file, 'r') as file:
        data = json.load(file)
        file_patterns = data.get('OUTPUT_FILE_PATTERNS', [])  # Provide default empty list if key doesn't exist
        print(f"File patterns for reduce: {file_patterns}")

for file_pattern in file_patterns:
    output_file_paths = get_outputfile_paths(sim_path, file_pattern)
    if len(output_file_paths) > 0:
        if output_file_paths[0].suffix == ".csv":
            merge_csv(output_file_paths=output_file_paths, output_path=output_dir)
        if (output_file_paths[0].suffix == ".phsp") or (output_file_paths[0].suffix == ".header"):
            phsp_files = [f for f in output_file_paths if f.suffix == ".phsp"]
            merge_phsp_files(phsp_files, output_dir)
            header_files = [f for f in output_file_paths if f.suffix == ".header"]
            merge_header_files(header_files, output_dir)

