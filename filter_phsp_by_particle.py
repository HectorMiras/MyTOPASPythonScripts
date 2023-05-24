import sys
import os
import os.path
import numpy as np
from phsp_manager import filter_phsp_by_particle

folder = '../tests/phsp_split'
name = 'nucleus_cylinder_PHSP_1mgml_0AuNP'
particle = 'electrons'
filter_phsp_by_particle(folder, name, particle)