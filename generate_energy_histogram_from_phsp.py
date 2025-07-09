import sys
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from phsp_manager import generate_weighted_energy_histogram

# folder = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/smallcellresults/topas-azure-med0p05-cells0p05/results'
folder = './../TOPAS_CellsNPs/work/topas-med1-smallcell1/results'
filename = 'PhaseSpace_NP' # without extension
min_energy = 0
max_energy = 0.036
bins = 100
particle = 'gammas' # 'electrons', 'gammas'
output_txt = True
show = True
save_path = True
generate_weighted_energy_histogram(folder, filename, bins, min_energy, max_energy, particle, show=show, save_path=save_path, output_txt=output_txt, log=False)

