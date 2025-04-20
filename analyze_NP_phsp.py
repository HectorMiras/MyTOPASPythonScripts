from phsp_manager import generate_weighted_energy_histogram, np_stats

folder = '/home/radiofisica/hector/mytopassimulations/TOPAS_CellsNPs/work/test-med1-smallcell1/results'
filename = 'PhaseSpace_NP' # without extension


min_energy = 0
max_energy = 0.025
bins = 1000
particle = 'electrons' # 'electrons', 'gammas'
output_txt = True
show = False
save_path = True
generate_weighted_energy_histogram(folder, filename, bins, min_energy, max_energy, particle, show, save_path, output_txt)

min_energy = 0
max_energy = 0.020
bins = 1000
particle = 'gammas' # 'electrons', 'gammas'
#generate_weighted_energy_histogram(folder, filename, bins, min_energy, max_energy, particle, show, save_path, output_txt)

np_stats(folder,filename)
