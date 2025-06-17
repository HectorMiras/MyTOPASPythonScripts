#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/26/23 3:48 PM

Script to show how to use the repair module to simulate variable dose rate functions

@author: alejandrobertolet
"""

import os, sys

# Add both the parent directory and the ChronoDNARepair directory to the Python path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
chrono_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../ChronoDNARepair'))
sys.path.insert(0, parent_dir)  # Add parent dir first
sys.path.insert(0, chrono_dir)  # Add ChronoDNARepair dir second

# Import our custom simulator that handles cell/run directory structure
from ChronoDNARepair.repair.custom_simulator import CustomSimulator

##############################
# SETUP OF REPAIR SIMULATION #
##############################

# Time options is a list with the initial time, final time and number of steps (or a list of custom time points as 4th arg)
# Times need to be given in seconds
initialTime = 0
finalTime = 48 * 3600 
nSteps = 48
timeOptions = [initialTime, finalTime, nSteps]

# Nucleus size in microns
nucleusMaxRadius = 4.65

# Models used for the simulation
diffusionModel = 'free'
dsbModel = 'standard'
ssbModel = 'standard'
bdModel = 'standard'

# Number of identical cells simulated
nCells = 40

# Dose rate function
doseratefunction = 'exponential'
initialDoseRate = 0.13803/3600 #  Gy/h
halfLife = (59.39*24) * 3600 #  59.39 days in seconds
# irradiationTime indicates the maximum length of exposure, it has to be provided to the simulator, otherwise it
# interprets instantaneous exposures
irradiationTime = 23 * 3600 # 23 hours

###############
# READ DAMAGE #
###############

# Set base path for SDD files with damage induced and dose to be loaded
# This path should contain cell# directories, which in turn contain run# directories with damage files
#damagepath = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/only_results_CellColony-med1-cell1/'
damagepath = '/home/radiofisica/hector/mytopassimulations/TOPAS_CellsNPs/work/CellColony-med0-cell0/'
maximumDose = -1 # Gy # This is a limit that is not used if the accumulated dose does not reach it

########################
# INITIALIZE SIMULATOR #
########################
# Simulator uses the previously defined options to initialize the simulation
# doseratefunctionargs is a list with the arguments of the dose rate function
# irradiationTime indicates the maximum length of exposure, it has to be provided to the simulator, otherwise it
# interprets instantaneous exposures

# Default cell parameters for the simulation 
cellpars = {
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

# Get and print DSB remaining output
dsbOutput = sim.avgRemainingDSBOverTime
times = dsbOutput.times
avgDSBremaining = dsbOutput.avgyvalues
for t in range(len(times)):
    print('Time: ', times[t], 'h, Fraction of DSB remaining: ', avgDSBremaining[t])

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

print("--------------------------------")

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
    report.append(f"Time: {times[t]:.1f}h, Fraction of DSB remaining: {avgDSBremaining[t]:.4f}")

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
report_path = os.path.join(damagepath, 'chronorepair_report.txt')
with open(report_path, 'w') as f:
    f.write('\n'.join(report))

# Print report to screen
print("\nFull Report:")
print("============")
for line in report:
    print(line)
print(f"\nReport saved to: {report_path}")
