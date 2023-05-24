import sys
import os
import os.path
import numpy as np
from phsp_manager import read_phsp_to_dataframe, plot_phsp_particle_positions
import matplotlib.pyplot as plt

folder = '/home/hector/mytopassimulations/MGHsimulations/tests/phsp_visualization'
name = 'nucleus_PHSP_1mgml_400000AGuIX_electrons'

data = read_phsp_to_dataframe(os.path.join(folder, name))
plot_phsp_particle_positions(data, os.path.join(folder, name+'.png'))

plt.show()