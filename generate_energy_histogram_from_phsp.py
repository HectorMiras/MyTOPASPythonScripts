import sys
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from phsp_manager import generate_weighted_energy_histogram

folder = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/smallcellresults/topas-azure-med1-cells1-phspnp/nodes_output/results'
filename = 'PhaseSpace_NP' # without extension
min_energy = 0
max_energy = 0.035
bins = 50
particle = 'gammas' # 'electrons', 'gammas'
output_txt = True
show = True
save_path = True
generate_weighted_energy_histogram(folder, filename, bins, min_energy, max_energy, particle, show, save_path, output_txt)

