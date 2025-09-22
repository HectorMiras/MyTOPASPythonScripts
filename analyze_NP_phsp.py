from phsp_manager import generate_weighted_energy_histogram, np_stats

folder = '/home/radiofisica/hector/mytopassimulations/TOPAS_CellsNPs/work/CellColony-med5-cell5/results'
filename = 'PhaseSpace_NP_all_cells' # without extension


min_energy = 0
max_energy = 0.035
bins = 100
particle = 'electrons' # 'electrons', 'gammas'
output_txt = True
show = True
save_path = True
log=True
relative=True
generate_weighted_energy_histogram(folder, filename, bins, min_energy, max_energy, particle, show, save_path, output_txt, log=log, relative=relative)

min_energy = 0
max_energy = 0.015
bins = 1000
particle = 'gammas' # 'electrons', 'gammas'
log=Truerelative=True
generate_weighted_energy_histogram(folder, filename, bins, min_energy, max_energy, particle, show, save_path, output_txt, log=log, relative=relative)

np_stats(folder,filename)
