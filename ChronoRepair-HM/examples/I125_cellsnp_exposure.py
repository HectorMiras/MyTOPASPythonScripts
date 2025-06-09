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
from custom_simulator import CustomSimulator

##############################
# SETUP OF REPAIR SIMULATION #
##############################

# Time options is a list with the initial time, final time and number of steps (or a list of custom time points as 4th arg)
# Times need to be given in seconds
initialTime = 0
finalTime = 23 * 3600 # 23 hours
nSteps = 23
timeOptions = [initialTime, finalTime, nSteps]

# Nucleus size in microns
nucleusMaxRadius = 4.65

# Models used for the simulation
diffusionModel = 'free'
dsbModel = 'standard'
ssbModel = 'standard'
bdModel = 'standard'

# Number of identical cells simulated
nCells = 10

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
damagepath = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/only_results_CellColony-med1-cell1/'
maximumDose = -1 # Gy # This is a limit that is not used if the accumulated dose does not reach it

########################
# INITIALIZE SIMULATOR #
########################
# Simulator uses the previously defined options to initialize the simulation
# doseratefunctionargs is a list with the arguments of the dose rate function
# irradiationTime indicates the maximum length of exposure, it has to be provided to the simulator, otherwise it
# interprets instantaneous exposures
sim = CustomSimulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel, bdmodel=bdModel,
                nucleusMaxRadius=nucleusMaxRadius, doseratefunction=doseratefunction, doseratefunctionargs=[initialDoseRate, halfLife],
                irradiationTime=irradiationTime)
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
print("--------------------------------")