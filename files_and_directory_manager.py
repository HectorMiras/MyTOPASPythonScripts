import os
from pathlib import Path
import glob
import re

def remove_part_suffix(filename):
    # Use regex to remove '_part#' from the filename
    new_filename = re.sub(r'_part\d+', '', filename)
    return new_filename


def get_outputfile_paths(parent_dir, output_pattern):
    # List to store the full path names of all output files
    output_file_paths = []

    # Loop through the subdirectories in the parent directory
    for subdir in os.listdir(parent_dir):
        if subdir.startswith('run'):
            # Find all output files in the current run folder that start with output_name
            run_folder_path = os.path.join(parent_dir, subdir)
            output_files = glob.glob(f"{run_folder_path}/{output_pattern.split('/')[-1]}")
            # Add the full path of each output file to the list
            for file_path in output_files:
                output_file_paths.append(Path(file_path).absolute())

    return output_file_paths