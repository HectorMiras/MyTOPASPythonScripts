import sys
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from phsp_manager import generate_weighted_energy_histogram

folder = './work/run1'
filename = 'nucleus_PHSP_1mgml_0AGuIX_electrons' # without extension
min_energy = 0
max_energy = 0.22
bins = 220
particle = 'electrons' # 'electrons', 'gammas'
output_txt = True
show = False
save_path = True
generate_weighted_energy_histogram(folder, filename, bins, min_energy, max_energy, particle, show, save_path, output_txt)

