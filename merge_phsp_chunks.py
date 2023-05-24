import sys
import os
import os.path
import numpy as np
from phsp_manager import merge_phsp_files, merge_header_files
from files_and_directory_manager import get_outputfile_paths


folder_path = '../tests/phsp_split/work'
name = "nucleus_cylinder_PHSP_1mgml_0AuNP"

extension = 'phsp'
output_file_paths = get_outputfile_paths(folder_path,name,extension)
output_path = '../tests/phsp_split/work'
merge_phsp_files(output_file_paths, f'{output_path}/{name}.{extension}')

extension = 'header'
output_file_paths = get_outputfile_paths(folder_path,name,extension)
output_path = '../tests/phsp_split/work'
merge_header_files(output_file_paths, f'{output_path}/{name}.{extension}')
