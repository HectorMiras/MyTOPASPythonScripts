#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/26/23 3:48 PM

Script to show how to use the repair module to simulate variable dose rate functions

@author: alejandrobertolet
"""

import os, sys
import numpy as np


# Add both the parent directory and the ChronoDNARepair directory to the Python path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
chrono_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../ChronoDNARepair'))
sys.path.insert(0, parent_dir)  # Add parent dir first
sys.path.insert(0, chrono_dir)  # Add ChronoDNARepair dir second

# Import our custom simulator that handles cell/run directory structure
from ChronoDNARepair.repair.custom_simulator import CustomSimulator
from my_chronorepair_plots import plot_timeofdeath_histogram

# Number of identical cells simulated
nCells = 80

###############
# READ DAMAGE #
###############

# Set base path for SDD files with damage induced and dose to be loaded
# This path should contain cell# directories, which in turn contain run# directories with damage files
#damagepath = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/only_results_CellColony-med1-cell1/'
damagepath = '/home/radiofisica/hector/mytopassimulations/TOPAS_CellsNPs/work/CellColony-med1-cell1/'
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


sim = CustomSimulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel, bdmodel=bdModel,
                nucleusMaxRadius=nucleusMaxRadius, doseratefunction=doseratefunction, doseratefunctionargs=[initialDoseRate, halfLife],
                irradiationTime=irradiationTime, cellparams=cellpars)
# Note: We don't need to call ReadDamage explicitly for a single cell
# because our custom Run method will handle reading damage from each cell directory

##################
# RUN SIMULATION #
##################
# This will use the custom Run method that reads damage from separate cell directories
# rereadDamageForNewRuns=True causes the simulator to look for different cell directories
# for each simulated cell, helping to represent cell-to-cell variation
sim.Run(nCells, rereadDamageForNewRuns=True, basepath=damagepath, maxDose=maximumDose, verbose=1, getVideo=False)

# Ensure chronorepair output directory exists
chronorepair_dir = os.path.join(damagepath, 'chronorepair')
os.makedirs(chronorepair_dir, exist_ok=True)

# Get and print DSB remaining output
dsbOutput = sim.avgRemainingDSBOverTime
times = dsbOutput.times
avgDSBremaining = dsbOutput.avgyvalues
varDSBremaining = dsbOutput.varyvalues
for t in range(len(times)):
    print('Time: ', times[t], 'h, Fraction of DSB remaining: ', avgDSBremaining[t], ' +/- ', np.sqrt(varDSBremaining[t]))

# Calculate and report cell survival fraction
celloutput = sim.celloutput
surviving_cells = sum(1 for cell in celloutput.celllist if cell.Surviving)
total_cells = len(celloutput.celllist)
survival_fraction = surviving_cells / total_cells if total_cells > 0 else 0

print("\n----- Cell Survival Results -----")
print(f"Surviving cells: {surviving_cells}/{total_cells}")
print(f"Survival fraction: {survival_fraction:.4f} ({surviving_cells/total_cells*100:.1f}%)")

# Calculate and report death causes
dead_cells = [cell for cell in celloutput.celllist if not cell.Surviving]
death_causes = {}
for cell in dead_cells:
    cause = cell._causeofdeath if hasattr(cell, '_causeofdeath') and cell._causeofdeath else "Unknown"
    death_causes[cause] = death_causes.get(cause, 0) + 1

if dead_cells:
    print("\n----- Causes of Cell Death -----")
    for cause, count in death_causes.items():
        percentage = (count / len(dead_cells)) * 100
        print(f"{cause}: {count} cells ({percentage:.1f}% of dead cells)")

    # --- CSV output for dead cells ---
    import csv
    csv_path = os.path.join(chronorepair_dir, 'chronorepair_cell_deaths.csv')
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Cell Index', 'Time of Death (h)', 'Type of Death'])
        for idx, cell in enumerate(celloutput.celllist, start=1):
            if not cell.Surviving:
                cause = getattr(cell, '_causeofdeath', None) or 'Unknown'
                # Use the output property, which is now set correctly
                time_of_death = getattr(cell, 'TimeOfDeath', None)
                if time_of_death is not None and float(time_of_death) >= 0:
                    try:
                        time_of_death_h = float(time_of_death) / 3600
                    except Exception:
                        time_of_death_h = 'N/A'
                else:
                    time_of_death_h = 'N/A'
                writer.writerow([idx, time_of_death_h, cause])
    print(f"\nCell death CSV saved to: {csv_path}")

print("--------------------------------")

plot_timeofdeath_histogram(csv_path, nbin=5, causes=["MitoticCatastrophe", "Apoptosis", "Necrosis"])

# Generate report
report = [
    "ChronoDNARepair Simulation Report",
    "==============================",
    "",
    "Simulation Parameters:",
    "--------------------",
    f"Time Range: {initialTime/3600:.1f}h to {finalTime/3600:.1f}h",
    f"Number of Time Steps: {nSteps}",
    f"Nucleus Max Radius: {nucleusMaxRadius} microns",
    f"Number of Cells: {nCells}",
    "",
    "Models:",
    "-------",
    f"Diffusion Model: {diffusionModel}",
    f"DSB Model: {dsbModel}",
    f"SSB Model: {ssbModel}",
    f"BD Model: {bdModel}",
    "",
    "Dose Rate Parameters:",
    "-------------------",
    f"Function Type: {doseratefunction}",
    f"Initial Dose Rate: {initialDoseRate*3600:.4f} Gy/h",
    f"Half Life: {halfLife/3600/24:.2f} days",
    f"Irradiation Time: {irradiationTime/3600:.1f}h",
    "",
    "DSB Repair Results:",
    "-----------------"
]

# Add DSB remaining data
for t in range(len(times)):
    report.append(f"Time: {times[t]:.1f}h, Fraction of DSB remaining: {avgDSBremaining[t]:.4f} +/- {np.sqrt(varDSBremaining[t]):.4f}")

report.extend([
    "",
    "Cell Survival Results:",
    "--------------------",
    f"Surviving cells: {surviving_cells}/{total_cells}",
    f"Survival fraction: {survival_fraction:.4f} ({surviving_cells/total_cells*100:.1f}%)",
])

# Add death causes if there are any dead cells
if dead_cells:
    report.extend([
        "",
        "Causes of Cell Death:",
        "------------------"
    ])
    for cause, count in death_causes.items():
        percentage = (count / len(dead_cells)) * 100
        report.append(f"{cause}: {count} cells ({percentage:.1f}% of dead cells)")

# Save report to file
report_path = os.path.join(chronorepair_dir, f'chronorepair_report_custom_params{nCells}_cells.txt')
with open(report_path, 'w') as f:
    f.write('\n'.join(report))

# Print report to screen
print("\nFull Report:")
print("============")
for line in report:
    print(line)
print(f"\nReport saved to: {report_path}")

