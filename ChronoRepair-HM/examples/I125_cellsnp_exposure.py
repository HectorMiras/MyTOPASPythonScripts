#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/26/23 3:48 PM

Script to show how to use the repair module to simulate variable dose rate functions

@author: alejandrobertolet
"""

import os, sys
import numpy as np
import csv


# Add both the parent directory and the ChronoDNARepair directory to the Python path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
chrono_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../ChronoDNARepair'))
sys.path.insert(0, parent_dir)  # Add parent dir first
sys.path.insert(0, chrono_dir)  # Add ChronoDNARepair dir second

# Import our custom simulator that handles cell/run directory structure
from ChronoDNARepair.repair.custom_simulator import CustomSimulator
from my_chronorepair_plots import plot_timeofdeath_histogram, plot_DSBremaining

# Number of identical cells simulated
nCells = 80

# Number of simulation repeats for statistics
Nrepeats = 10  # Set as desired

###############
# READ DAMAGE #
###############

# Set base path for SDD files with damage induced and dose to be loaded
# This path should contain cell# directories, which in turn contain run# directories with damage files
#damagepath = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/only_results_CellColony-med1-cell1/'
damagepath = '/home/radiofisica/hector/mytopassimulations/TOPAS_CellsNPs/work/CellColony-med0-cell0/'
maximumDose = -1 # Gy # This is a limit that is not used if the accumulated dose does not reach it

def read_total_avg_dose_per_cell(damagepath, nCells):
    """
    Reads the total delivered dose for each cell (by summing all run# DNADamage.phsp files)
    and returns the average dose per cell.
    """
    total_doses = []
    for cell_idx in range(1, nCells+1):
        cell_dir = os.path.join(damagepath, f'cell{cell_idx}')
        if not os.path.isdir(cell_dir):
            print(f"Warning: {cell_dir} does not exist, skipping.")
            continue
        cell_dose = 0.0
        for run_name in os.listdir(cell_dir):
            run_dir = os.path.join(cell_dir, run_name)
            if not (run_name.startswith('run') and os.path.isdir(run_dir)):
                continue
            phsp_path = os.path.join(run_dir, 'DNADamage.phsp')
            if not os.path.isfile(phsp_path):
                continue
            with open(phsp_path, 'r') as f:
                lines = f.readlines()
                if not lines:
                    continue
                # The last line's second column is the cumulative dose for this run
                try:
                    last_line = lines[-1]
                    split = last_line.split()
                    run_dose = float(split[1])
                    cell_dose += run_dose
                except Exception as e:
                    print(f"Error reading dose from {phsp_path}: {e}")
        total_doses.append(cell_dose)
    if total_doses:
        avg_dose = sum(total_doses) / len(total_doses)
    else:
        avg_dose = 0.0
    print(f"Average dose per cell: {avg_dose} Gy")
    return avg_dose

D_total_from_phsps = read_total_avg_dose_per_cell(damagepath, nCells)

##############################
# SETUP OF REPAIR SIMULATION #
##############################

# Time options is a list with the initial time, final time and number of steps (or a list of custom time points as 4th arg)
# Times need to be given in seconds
initialTime = 0
# finalTime = 168 * 3600 
# nSteps = 48
# timeOptions = [initialTime, finalTime, nSteps]

time_custom=[]
time_custom +=[t*3600 for t in range(0,49,1)] # Every hour for the first 48 hours
time_custom +=[t*3600 for t in range(51,121,3)] # Every 3 hours from 51 to 121 hours
time_custom +=[t*3600 for t in range(132,241,12)] # Every 12 hours from 132 to 241 hours  
finalTime = time_custom[-1]
nSteps = len(time_custom)
timeOptions = [initialTime, 240*3600, len(time_custom), time_custom]


# Nucleus size in microns
nucleusMaxRadius = 4.65

# Models used for the simulation
diffusionModel = 'free'
dsbModel = 'standard'
ssbModel = 'standard'
bdModel = 'standard'


# Dose rate function
doseratefunction = 'exponential' # 'uniform', 'exponential', 'linear'
# Dose_without_NPs = 3.55 # Gy
# dose_enhancement = 1.16 # Dose enhancement factor due to NPs 1mg/ml=1.16
# D_total = Dose_without_NPs * dose_enhancement # Total dose to be delivered, Gy
D_total = D_total_from_phsps

halfLife = (59.39*24) * 3600 #  59.39 days in seconds
lam = np.log(2) / halfLife
# irradiationTime indicates the maximum length of exposure, it has to be provided to the simulator, otherwise it
# interprets instantaneous exposures
irradiationTime = 24 * 3600 # 24 hours

initialDoseRate = D_total * lam / (1 - np.exp(-lam*irradiationTime)) #  Gy/s
print("Initial dose rate (Gy/s): ", initialDoseRate)




########################
# INITIALIZE SIMULATOR #
########################
# Simulator uses the previously defined options to initialize the simulation
# doseratefunctionargs is a list with the arguments of the dose rate function
# irradiationTime indicates the maximum length of exposure, it has to be provided to the simulator, otherwise it
# interprets instantaneous exposures

# Default cell parameters for the simulation 
cellpars_default = {
    'cycling': True,
    'initialStage': 'G1',
    'mu_G1': 10.0, 'sigma_G1': 0.1, 'unit_G1': 'h',
    'mu_S': 7.0, 'sigma_S': 0.1, 'unit_S': 'h',
    'mu_G2': 5.0, 'sigma_G2': 0.1, 'unit_G2': 'h',
    'mu_M': 1.5, 'sigma_M': 0.1, 'unit_M': 'h',
    'propHeterochromatin_G1': 0.1,
    'propHeterochromatin_S': 0.25,
    'propHeterochromatin_G2': 0.5,
    'propHeterochromatin_M': 0.9,
    'damageThresholdForNecrosisAtG1S': 20,
    'damageThresholdForApoptosisAtG1S': 14,
    'damageThresholdForArrestAtG1S': 8,
    'damageThresholdForNecrosisAtG2M': 16,
    'damageThresholdForApoptosisAtG2M': 8,
    'damageThresholdForArrestAtG2M': 5,
    'misrepairThresholdForMitoticCatastrophe': 5,
    'mu_necrosis': 2.0, 'sigma_necrosis': 0.75, 'unit_necrosis': 'h',
    'mu_apoptosis': 3.0, 'sigma_apoptosis': 0.5, 'unit_apoptosis': 'h',
    'mu_mitoticcat': 10.0, 'sigma_mitoticcat': 0.75, 'unit_mitoticcat': 'h',
    'oxygen_concentration': 1.0
}

# custom params to fit observed survival
A_th = 1.2
MC_th = 1.2
# r_rep = 1.15

cellpars = {
    'cycling': True,
    'initialStage': 'G1',
    'mu_G1': 27.7, 'sigma_G1': 0.1, 'unit_G1': 'h',
    'mu_S': 19.3, 'sigma_S': 0.1, 'unit_S': 'h',
    'mu_G2': 13.8, 'sigma_G2': 0.1, 'unit_G2': 'h',
    'mu_M': 4.15, 'sigma_M': 0.1, 'unit_M': 'h',
    'propHeterochromatin_G1': 0.1,
    'propHeterochromatin_S': 0.25,
    'propHeterochromatin_G2': 0.5,
    'propHeterochromatin_M': 0.9,
    'damageThresholdForNecrosisAtG1S': 20,
    'damageThresholdForApoptosisAtG1S': 14*A_th,
    'damageThresholdForArrestAtG1S': 8,
    'damageThresholdForNecrosisAtG2M': 16,
    'damageThresholdForApoptosisAtG2M': 8*A_th,
    'damageThresholdForArrestAtG2M': 5,
    'misrepairThresholdForMitoticCatastrophe': 5*MC_th,
    'mu_necrosis': 2.0, 'sigma_necrosis': 0.75, 'unit_necrosis': 'h',
    'mu_apoptosis': 3.0, 'sigma_apoptosis': 0.5, 'unit_apoptosis': 'h',
    'mu_mitoticcat': 10.0, 'sigma_mitoticcat': 0.75, 'unit_mitoticcat': 'h',
    'oxygen_concentration': 1.0
}

# Ensure chronorepair output directory exists
chronorepair_dir = os.path.join(damagepath, 'chronorepair')
os.makedirs(chronorepair_dir, exist_ok=True)

# For statistics
all_total_deaths = []
all_death_cause_counts = []
all_survival_fractions = []
all_total_cells = []
all_avgDSBremaining = []
all_varDSBremaining = []
all_times = []
all_time_of_death = []
all_dead_cell_info = []  # List of lists: for each repeat, a list of (cell_idx, time_of_death_h, cause)
all_DSBremaining_per_cell = []  # List of arrays: for each repeat, shape (nCells, nTimePoints)

for rep in range(Nrepeats):
    print(f"\n=== Simulation Repeat {rep+1}/{Nrepeats} ===")
    sim = CustomSimulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel, bdmodel=bdModel,
                    nucleusMaxRadius=nucleusMaxRadius, doseratefunction=doseratefunction, doseratefunctionargs=[initialDoseRate, halfLife],
                    irradiationTime=irradiationTime, cellparams=cellpars)
    sim.Run(nCells, rereadDamageForNewRuns=True, basepath=damagepath, maxDose=maximumDose, verbose=0, getVideo=False)
    celloutput = sim.celloutput
    total_cells = len(celloutput.celllist)
    dead_cells = [cell for cell in celloutput.celllist if not cell.Surviving]
    total_deaths = len(dead_cells)
    survival_fraction = (total_cells - total_deaths) / total_cells if total_cells > 0 else 0
    # Count by cause
    death_causes = {}
    rep_dead_cell_info = []
    for idx, cell in enumerate(celloutput.celllist, start=1):
        if not cell.Surviving:
            cause = getattr(cell, '_causeofdeath', None) or 'Unknown'
            time_of_death = getattr(cell, 'TimeOfDeath', None)
            if time_of_death is not None and float(time_of_death) >= 0:
                try:
                    time_of_death_h = float(time_of_death) / 3600
                except Exception:
                    time_of_death_h = None
            else:
                time_of_death_h = None
            if time_of_death_h is not None:
                rep_dead_cell_info.append((idx, time_of_death_h, cause))
            death_causes[cause] = death_causes.get(cause, 0) + 1
    all_dead_cell_info.append(rep_dead_cell_info)
    all_total_deaths.append(total_deaths)
    all_death_cause_counts.append(death_causes)
    all_survival_fractions.append(survival_fraction)
    all_total_cells.append(total_cells)
    # DSB remaining output
    dsbOutput = sim.avgRemainingDSBOverTime
    all_times.append(np.array(dsbOutput.times))
    all_avgDSBremaining.append(np.array(dsbOutput.avgyvalues))
    all_varDSBremaining.append(np.array(dsbOutput.varyvalues))
    # Collect per-cell DSB remaining for this repeat (if available)
    if hasattr(dsbOutput, 'runlist') and len(dsbOutput.runlist) == nCells:
        # Each runlist[i] is a TimeCurveForSingleRun for cell i
        per_cell = np.array([tc.yvalues for tc in dsbOutput.runlist])  # shape (nCells, nTimePoints)
        all_DSBremaining_per_cell.append(per_cell)
    # Time of death for all dead cells in this repeat
    rep_time_of_death = [info[1] for info in rep_dead_cell_info]
    all_time_of_death.append(rep_time_of_death)

# Aggregate DSB remaining (mean and std over all cells in all repeats)
if all_DSBremaining_per_cell:
    all_DSB_cells = np.concatenate(all_DSBremaining_per_cell, axis=0)  # shape (Nrepeats*nCells, nTimePoints)
    mean_DSB_cells = np.mean(all_DSB_cells, axis=0)
    std_DSB_cells = np.std(all_DSB_cells, axis=0)
    # Also calculate std of the mean (std of per-repeat means)
    per_repeat_means = np.array([np.mean(per_cell, axis=0) for per_cell in all_DSBremaining_per_cell])  # shape (Nrepeats, nTimePoints)
    std_mean_DSB_cells = np.std(per_repeat_means, axis=0)
    ref_times = all_times[0]
else:
    mean_DSB_cells = np.zeros_like(all_avgDSBremaining[0])
    std_DSB_cells = np.zeros_like(all_avgDSBremaining[0])
    std_mean_DSB_cells = np.zeros_like(all_avgDSBremaining[0])
    ref_times = all_times[0]

# Save DSB remaining mean/std as CSV (cell-to-cell variability and std of mean)
dsb_csv_path = os.path.join(chronorepair_dir, f'chronorepair_DSBremaining_{nCells}cells_{Nrepeats}repeats.csv')
with open(dsb_csv_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Time (h)', 'Mean Fraction DSB Remaining', 'Std Fraction DSB Remaining', 'Std_mean Fraction DSB Remaining'])
    for t, mean_val, std_val, std_mean_val in zip(ref_times, mean_DSB_cells, std_DSB_cells, std_mean_DSB_cells):
        writer.writerow([t, mean_val, std_val, std_mean_val])
print(f"DSB remaining mean/std CSV saved to: {dsb_csv_path}")

# Aggregate time of death for all repeats (flattened)
flat_time_of_death = [t for rep_list in all_time_of_death for t in rep_list]

# Aggregate death causes
from collections import Counter, defaultdict

def aggregate_death_cause_stats(all_death_cause_counts):
    all_causes = set()
    for d in all_death_cause_counts:
        all_causes.update(d.keys())
    stats = {}
    for cause in all_causes:
        vals = [d.get(cause, 0) for d in all_death_cause_counts]
        stats[cause] = {
            'mean': np.mean(vals),
            'std': np.std(vals),
            'all': vals
        }
    return stats

death_cause_stats = aggregate_death_cause_stats(all_death_cause_counts)
mean_total_deaths = np.mean(all_total_deaths)
std_total_deaths = np.std(all_total_deaths)
mean_survival_fraction = np.mean(all_survival_fractions)
std_survival_fraction = np.std(all_survival_fractions)

# Reporting
print("\n===== Statistical Summary over", Nrepeats, "repeats =====")
print(f"Mean total cell deaths: {mean_total_deaths:.2f} ± {std_total_deaths:.2f}")
print(f"Mean survival fraction: {mean_survival_fraction:.4f} ± {std_survival_fraction:.4f}")
print("Deaths by cause (mean ± std):")
for cause, stats in death_cause_stats.items():
    print(f"  {cause}: {stats['mean']:.2f} ± {stats['std']:.2f}")

# Save reporting section to text file
stats_txt_path = os.path.join(chronorepair_dir, 'chronorepair_stats_summary.txt')
with open(stats_txt_path, 'w') as f:
    f.write(f"===== Statistical Summary over {Nrepeats} repeats =====\n")
    f.write(f"Mean total cell deaths: {mean_total_deaths:.2f} ± {std_total_deaths:.2f}\n")
    f.write(f"Mean survival fraction: {mean_survival_fraction:.4f} ± {std_survival_fraction:.4f}\n")
    f.write("Deaths by cause (mean ± std):\n")
    for cause, stats in death_cause_stats.items():
        f.write(f"  {cause}: {stats['mean']:.2f} ± {stats['std']:.2f}\n")
print(f"Statistical summary saved to: {stats_txt_path}")

# Prepare all_causes for summary table
all_causes = set()
for d in all_death_cause_counts:
    all_causes.update(d.keys())
all_causes = sorted(all_causes)
# Save summary CSV in the requested per-iteration table format (no printing)
summary_csv_path = os.path.join(chronorepair_dir, f'chronorepair_summary_{nCells}cells_{Nrepeats}repeats.csv')
with open(summary_csv_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    # Header
    header = ['Death type'] + [f'Iteration {i+1}' for i in range(Nrepeats)]
    writer.writerow(header)
    # Total deaths row
    writer.writerow(['Total'] + all_total_deaths)
    # Per-cause rows
    for cause in all_causes:
        row = [cause]
        for d in all_death_cause_counts:
            row.append(d.get(cause, 0))
        writer.writerow(row)

# Save all time of death values for all repeats in the same structure as the provided CSV
all_timeofdeath_csv = os.path.join(chronorepair_dir, f'chronorepair_timeofdeath_{nCells}cells_{Nrepeats}repeats.csv')
with open(all_timeofdeath_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Cell Index', 'Time of Death (h)', 'Type of Death'])
    for rep_dead_cell_info in all_dead_cell_info:
        for idx, time_of_death_h, cause in rep_dead_cell_info:
            writer.writerow([idx, time_of_death_h, cause])
print(f"All time of death values saved to: {all_timeofdeath_csv}")

# Optionally, plot DSB remaining
plot_DSBremaining(dsb_csv_path, std_type='std')
# Optionally, plot histogram of all time of death values
plot_timeofdeath_histogram(all_timeofdeath_csv, nbin=6)

