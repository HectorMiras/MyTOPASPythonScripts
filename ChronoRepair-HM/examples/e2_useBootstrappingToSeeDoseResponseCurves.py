#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/26/23 1:50 PM

This script shows how to use the bootstrapping method to see the dose response curves with better statistics

@author: alejandrobertolet
"""

import os, sys
# Add both the parent directory and the ChronoDNARepair directory to the Python path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
chrono_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../ChronoDNARepair'))
sys.path.insert(0, parent_dir)  # Add parent dir first
sys.path.insert(0, chrono_dir)  # Add ChronoDNARepair dir second
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from ChronoDNARepair.induction.damage import DamageToDNA

#####################################
# READING THE DAMAGE FROM SDD FILES #
#####################################

# Set base path for SDD files with damage induced and dose to be loaded
damagepath = './damageFromTopas-nBio/xray-250keV/'
maxDose = 0.5 # Gy

# Read damage
nfiles = len(os.listdir(damagepath))

# Section to get what directories actually contains both dose and SDD. Disregard others!
listOfAvailableDirs = []
for j in range(nfiles):
    newpath = damagepath + str(j) + '/'
    if str(j) in os.listdir(damagepath):
        files = os.listdir(newpath)
        if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
            if os.path.getsize(newpath + 'DNADamage_sdd.txt') > 0:  # only those with actual data
                listOfAvailableDirs.append(j)

# Preparing arrays for multiple damage readouts
Dose = np.linspace(0, maxDose, 100)
DSB = np.zeros(Dose.shape)
nSites = np.zeros(Dose.shape)
SSB = np.zeros(Dose.shape)
BD = np.zeros(Dose.shape)
ntries = np.zeros(Dose.shape)

#################
# BOOTSTRAPPING #
#################
nboot = 5 # number of bootstraps
# Prepare arrays for bootstrapping
alldsbs = np.zeros([len(Dose), nboot])
allssbs = np.zeros([len(Dose), nboot])
allbds = np.zeros([len(Dose), nboot])
allnsites = np.zeros([len(Dose), nboot])
for i in range(nboot):
    print('Bootstrap number: ' + str(i + 1))
    # Initialize damage
    damage = DamageToDNA()
    # Get new random sequence of directories
    neworder = np.random.choice(listOfAvailableDirs, len(listOfAvailableDirs), replace=False)
    # Read SDD files
    for j, e in enumerate(neworder):
        path = damagepath + str(e) + '/'
        damage.readSDDAndDose(path)

    # Populate damage object with damage induced
    damage.populateDamages(getVideo=False, stopAtDose=maxDose, stopAtTime=0.0, recalculatePerEachTrack=True)
    # Get dose-response curve
    dose, dsb = damage.getDoseResponseCurve(q='DSB', plot=False)
    dose, ssb = damage.getDoseResponseCurve(q='SSB', plot=False)
    dose, bd = damage.getDoseResponseCurve(q='BD', plot=False)
    dose, nsites = damage.getDoseResponseCurve(q='nSites', plot=False)
    # Convert into numpy arrays
    dose = np.array(dose)
    dsb = np.array(dsb)
    ssb = np.array(ssb)
    bd = np.array(bd)
    nsites = np.array(nsites)
    # Add zeros
    dose = np.insert(dose, 0, 0)
    dsb = np.insert(dsb, 0, 0)
    ssb = np.insert(ssb, 0, 0)
    bd = np.insert(bd, 0, 0)
    nsites = np.insert(nsites, 0, 0)
    # Get interpolation of dose-response curve
    fdsb = interpolate.interp1d(dose, dsb, kind='linear')
    fssb = interpolate.interp1d(dose, ssb, kind='linear')
    fbd = interpolate.interp1d(dose, bd, kind='linear')
    fnsites = interpolate.interp1d(dose, nsites, kind='linear')
    # Get dose-response curve for the dose range
    for id, d in enumerate(Dose):
        if d >= 0 and d < np.max(dose):
            DSB[id] += fdsb(d)
            SSB[id] += fssb(d)
            BD[id] += fbd(d)
            nSites[id] += fnsites(d)
            alldsbs[id, i] = fdsb(d)
            allssbs[id, i] = fssb(d)
            allbds[id, i] = fbd(d)
            allnsites[id, i] = fnsites(d)
            ntries[id] += 1

# Get statistics
vardsb = np.zeros(Dose.shape)
varssb = np.zeros(Dose.shape)
varbd = np.zeros(Dose.shape)
varnsites = np.zeros(Dose.shape)
for j in range(alldsbs.shape[0]):
    vardsb[j] = np.var(alldsbs[j, :])
    varssb[j] = np.var(allssbs[j, :])
    varbd[j] = np.var(allbds[j, :])
    varnsites[j] = np.var(allnsites[j, :])
stddsb = np.sqrt(vardsb)
stdssb = np.sqrt(varssb)
stdbd = np.sqrt(varbd)
stdnsites = np.sqrt(varnsites)
# Get average
DSB /= ntries
SSB /= ntries
BD /= ntries
nSites /= ntries

################
# PLOT RESULTS #
################

fig, ax = plt.subplots(2, 2, figsize=(10, 10))
ax[0, 0].plot(Dose, DSB, 'k-', label='DSB')
ax[0, 0].fill_between(Dose, DSB - stddsb, DSB + stddsb, color='k', alpha=0.2)
ax[0, 0].set_xlabel('Dose (Gy)')
ax[0, 0].set_ylabel('Cumulative DSB')
ax[0, 0].grid()
ax[0, 0].set_xlim([0, maxDose * 1.01])
ax[0, 0].set_ylim([0, None])

ax[0, 1].plot(Dose, SSB, 'k-', label='SSB')
ax[0, 1].fill_between(Dose, SSB - stdssb, SSB + stdssb, color='k', alpha=0.2)
ax[0, 1].set_xlabel('Dose (Gy)')
ax[0, 1].set_ylabel('Cumulative SSB')
ax[0, 1].grid()
ax[0, 1].set_xlim([0, maxDose * 1.01])
ax[0, 1].set_ylim([0, None])

ax[1, 0].plot(Dose, BD, 'k-', label='BD')
ax[1, 0].fill_between(Dose, BD - stdbd, BD + stdbd, color='k', alpha=0.2)
ax[1, 0].set_xlabel('Dose (Gy)')
ax[1, 0].set_ylabel('Cumulative BD')
ax[1, 0].grid()
ax[1, 0].set_xlim([0, maxDose * 1.01])
ax[1, 0].set_ylim([0, None])

ax[1, 1].plot(Dose, nSites, 'k-', label='nSites')
ax[1, 1].fill_between(Dose, nSites - stdnsites, nSites + stdnsites, color='k', alpha=0.2)
ax[1, 1].set_xlabel('Dose (Gy)')
ax[1, 1].set_ylabel('Cumulative number of damage sites')
ax[1, 1].grid()
ax[1, 1].set_xlim([0, maxDose * 1.01])
ax[1, 1].set_ylim([0, None])

# Create temp directory if it doesn't exist
temp_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'temp')
os.makedirs(temp_dir, exist_ok=True)

plt.tight_layout()
plt.savefig(os.path.join(temp_dir, 'dose_response_curves.png'))
plt.close()
