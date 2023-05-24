import sys
import os
import os.path
import numpy as np
from phsp_manager import generate_header

folder = '../tests/phsp_split'
name = 'nucleus_cylinder_PHSP_1mgml_0AuNP'

generate_header(folder, name)