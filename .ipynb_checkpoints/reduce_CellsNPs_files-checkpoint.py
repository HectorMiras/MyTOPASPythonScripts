from topas_csv_files_manager import collect_np_number, merge_TOPAS_csv
from files_and_directory_manager import get_outputfile_paths
from phsp_manager import merge_phsp_files, merge_header_files
import json
import os
import sys
from pathlib import Path
import glob

sim_path = sys.argv[1]
sim_config_file = sys.argv[2]
#sim_path = "/home/hector/mytopassimulations/MGHsimulations/tests/merge_binned_outputs/"
#sim_config_file = "/home/hector/mytopassimulations/MGHsimulations/tests/merge_binned_outputs/simconfig.json"

print(f'Simulation path: {sim_path}')
print(f'Config file: {sim_config_file}')
#file_patterns = ["DoseToCell*", "DoseToNucleus*electrons.csv", "DoseToNucleus*gammas.csv", "nucleus_PHSP*"]

# If simconfig.json file exists in the directory, get the OUTPUT_FILE_PATTERNS from the file
if os.path.isfile(sim_config_file):
    # Open the file and read the JSON content
    with open(sim_config_file, 'r') as file:
        data = json.load(file)
        file_patterns = data.get('OUTPUT_FILE_PATTERNS', [])  # Provide default empty list if key doesn't exist
        print(f"File patterns for reduce: {file_patterns}")



for file_pattern in file_patterns:
    subdirectories = [d for d in Path(sim_path).iterdir() if d.is_dir() and d.name.startswith('run')]
    for rundir in subdirectories:
        output_file_paths = glob.glob(f"{rundir}/{file_pattern.split('/')[-1]}")
        # Filter the list to include only paths that contain "_part"
        filtered_paths = [Path(path) for path in output_file_paths if "_part" in str(path)]
        if len(filtered_paths) > 0:
            output_dir = str(Path(filtered_paths[0]).parent)
            if filtered_paths[0].suffix == ".csv":
                merge_TOPAS_csv(output_file_paths=filtered_paths, output_path=output_dir)
            if (filtered_paths[0].suffix == ".phsp") or (filtered_paths[0].suffix == ".header"):
                phsp_files = [f for f in filtered_paths if f.suffix == ".phsp"]
                merge_phsp_files(phsp_files, output_dir)
                header_files = [f for f in filtered_paths if f.suffix == ".header"]
                merge_header_files(header_files, output_dir)

for file_pattern in file_patterns:
    output_file_paths = get_outputfile_paths(sim_path, file_pattern)
    # Filter the list to include only paths not containing "_part"
    filtered_paths = [path for path in output_file_paths if "_part" not in str(path)]
    output_dir = os.path.join(sim_path, "results")
    if len(filtered_paths) > 0:
        if ('np_number' in str(filtered_paths[0])):
            collect_np_number(output_file_paths=filtered_paths, output_path=output_dir)
        if filtered_paths[0].suffix == ".csv":
            merge_TOPAS_csv(output_file_paths=filtered_paths, output_path=output_dir)
        if (filtered_paths[0].suffix == ".phsp") or (filtered_paths[0].suffix == ".header"):
            phsp_files = [f for f in filtered_paths if f.suffix == ".phsp"]
            merge_phsp_files(phsp_files, output_dir)
            header_files = [f for f in filtered_paths if f.suffix == ".header"]
            merge_header_files(header_files, output_dir)

