from topas_csv_files_manager import merge_combine_csv
from files_and_directory_manager import get_outputfile_paths
from phsp_manager import merge_phsp_files, merge_header_files
import json
import os
import sys

sim_path = sys.argv[1]
#sim_path = "/home/hector/mytopassimulations/MGHsimulations/tests/nodes_output"
output_dir = os.path.join(sim_path, "results")

file_patterns = ["DoseToCell*", "DoseToNucleus*electrons.csv", "DoseToNucleus*gammas.csv", "nucleus_PHSP*"]

for file_pattern in file_patterns:
    output_file_paths = get_outputfile_paths(sim_path, file_pattern)
    if len(output_file_paths) > 0:
        if output_file_paths[0].suffix == ".csv":
            merge_combine_csv(output_file_paths, output_dir)
        if (output_file_paths[0].suffix == ".phsp") or (output_file_paths[0].suffix == ".header"):
            phsp_files = [f for f in output_file_paths if f.suffix == ".phsp"]
            merge_phsp_files(phsp_files, output_dir)
            header_files = [f for f in output_file_paths if f.suffix == ".header"]
            merge_header_files(header_files, output_dir)

