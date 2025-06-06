#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/26/23 2:36 PM

This script shows how to simulate repair from an instantaneous exposure to a given dose of a given radiation

@author: alejandrobertolet
"""

import os, sys

# Add both the parent directory and the ChronoDNARepair directory to the Python path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
chrono_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../ChronoDNARepair'))
sys.path.insert(0, parent_dir)  # Add parent dir first
sys.path.insert(0, chrono_dir)  # Add ChronoDNARepair dir second
from ChronoDNARepair.repair.running import Simulator

##############################
# SETUP OF REPAIR SIMULATION #
##############################

# Time options is a list with the initial time, final time and number of steps (or a list of custom time points as 4th arg)
# Times need to be given in seconds
initialTime = 0
finalTime = 24 * 3600 # 24 hours
nSteps = 24
timeOptions = [initialTime, finalTime, nSteps]

# Nucleus size in microns
nucleusMaxRadius = 4.65

# Models used for the simulation
cellModel = 'standard'
diffusionModel = 'subdiffusion'
dsbModel = 'foci_cycle'
ssbModel = 'standard'
bdModel = 'standard'

# Number of identical cells simulated
nCells = 10

# No dose rate function is used
doseratefunction = None

###############
# READ DAMAGE #
###############

damagepath = './damageFromTopas-nBio/xray-250keV/'
instantaneousDose = 1.0 # Gy

########################
# INITIALIZE SIMULATOR #
########################
# Simulator uses the previously defined options to initialize the simulation
sim = Simulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel, bdmodel=bdModel,
                nucleusMaxRadius=nucleusMaxRadius, doseratefunction=doseratefunction, cellmodel=cellModel)
# Reads damage
sim.ReadDamage(damagepath, instantaneousDose)

##################
# RUN SIMULATION #
##################
# Options
# rereadDamageForNewRuns: If True, the damage is read again for each new cell. If False, the same initial damage is used for all cells
# basepath: Path to the damage files
# maxDose: Dose of the damage files
# version: Version of the damage files
# verbose: Verbosity level
# getVideo: If True, a video of the simulation is saved
sim.Run(nCells, rereadDamageForNewRuns=True, basepath=damagepath, maxDose=instantaneousDose, verbose=1, getVideo=False)

# Get and print output
dsbOutput = sim.avgRemainingDSBOverTime
times = dsbOutput.times
avgDSBremaining = dsbOutput.avgyvalues
for t in range(len(times)):
    print('Time: ', times[t], 'h, Fraction of DSB remaining: ', avgDSBremaining[t])