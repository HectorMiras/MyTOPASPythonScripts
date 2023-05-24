import sys
import os
import os.path
import numpy as np
from phsp_manager import split_phsp_file, generate_header

folder_path = '../tests/phsp_split/'
input_phsp_file = 'nucleus_cylinder_PHSP_1mgml_0AuNP'  # without extension
num_parts = 4  # specify the desired number of parts

split_phsp_file(folder_path, input_phsp_file, num_parts)
