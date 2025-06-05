#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/26/23 3:48 PM

Script to show how to use the repair module to simulate variable dose rate functions

@author: alejandrobertolet
"""

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
diffusionModel = 'free'
dsbModel = 'standard'
ssbModel = 'standard'
bdModel = 'standard'

# Number of identical cells simulated
nCells = 10

# Dose rate function
doseratefunction = 'exponential'
initialDoseRate = 0.1/3600 # 0.1 Gy/h
halfLife = 8 * 3600 # 8 hours
irradiationTime = 24 * 3600 # 24 hours

###############
# READ DAMAGE #
###############

damagepath = './damageFromTopas-nBio/xray-250keV/'
maximumDose = 10.0 # Gy # This is a limit that is not used if the accumulated dose does not reach it

########################
# INITIALIZE SIMULATOR #
########################
# Simulator uses the previously defined options to initialize the simulation
# doseratefunctionargs is a list with the arguments of the dose rate function
# irradiationTime indicates the maximum length of exposure, it has to be provided to the simulator, otherwise it
# interprets instantaneous exposures
sim = Simulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel, bdmodel=bdModel,
                nucleusMaxRadius=nucleusMaxRadius, doseratefunction=doseratefunction, doseratefunctionargs=[initialDoseRate, halfLife],
                irradiationTime=irradiationTime)
# Reads damage
sim.ReadDamage(damagepath, maximumDose)

##################
# RUN SIMULATION #
##################
sim.Run(nCells, rereadDamageForNewRuns=True, basepath=damagepath, maxDose=maximumDose, verbose=1, getVideo=False)

# Get and print output
dsbOutput = sim.avgRemainingDSBOverTime
times = dsbOutput.times
avgDSBremaining = dsbOutput.avgyvalues
for t in range(len(times)):
    print('Time: ', times[t], 'h, Fraction of DSB remaining: ', avgDSBremaining[t])