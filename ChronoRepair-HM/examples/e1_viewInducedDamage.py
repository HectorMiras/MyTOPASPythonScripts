#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/26/23 11:35 AM

This script shows how to use the DamageToDNA class to read SDD files and visualize the damage induced.

@author: alejandrobertolet
"""

import os, random
import sys

# Add both the parent directory and the ChronoDNARepair directory to the Python path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
chrono_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../ChronoDNARepair'))
sys.path.insert(0, parent_dir)  # Add parent dir first
sys.path.insert(0, chrono_dir)  # Add ChronoDNARepair dir second
from induction.damage import DamageToDNA  # Now Python will find it in ChronoDNARepair/induction/damage.py

#####################################
# READING THE DAMAGE FROM SDD FILES #
#####################################

# Set base path for SDD files with damage induced and dose to be loaded
#damagepath = './damageFromTopas-nBio/xray-250keV/'
damagepath = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/only_results_CellColony-med0-cell0/cell2/'
dose = -1  # Gy

# Define colormap for visualization
MyColorMap = 'viridis'

# Initialize damage object
damage = DamageToDNA()

# Helper function to get numbered directories
def get_numbered_dirs(path):
    """Get all numbered directories regardless of prefix"""
    numbered_dirs = []
    try:
        for d in os.listdir(path):
            full_path = os.path.join(path, d)
            if os.path.isdir(full_path):
                # Extract number from directory name using last digits
                number = ''.join(filter(str.isdigit, d))
                if number:  # Only add if there are digits in the name
                    numbered_dirs.append((int(number), d))
    except OSError:
        print(f"Error accessing directory: {path}")
        return []
    
    # Sort by the numeric value
    numbered_dirs.sort(key=lambda x: x[0])
    return [d[1] for d in numbered_dirs]  # Return only directory names

# Read SDD files
# Section to get directories that contain both dose and SDD files
listOfAvailableDirs = []
for dir_name in get_numbered_dirs(damagepath):
    newpath = os.path.join(damagepath, dir_name)
    files = os.listdir(newpath)
    if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
        if os.path.getsize(os.path.join(newpath, 'DNADamage_sdd.txt')) > 0:  # only those with actual data
            listOfAvailableDirs.append(dir_name)
neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
print(f"Number of directories with damage data: {len(neworder)}")
for i, e in enumerate(neworder):
    path = os.path.join(damagepath, e) + '/'
    #print(f"Reading damage from: {path}")
    # Read SDD file and dose
    damage.readSDDAndDose(path)

# Populate damage object with damage induced
# Options: getVideo to see how damage is populated track by track (then recalculatePerEachTrack has to be True);
# stopAtDose to get to certain dose; stopAtTime is only relevant when simulating repair or assigning times to damage
damage.populateDamages(getVideo=False, stopAtDose=dose, stopAtTime=0.0, recalculatePerEachTrack=False)

###############################
# METHODS TO VISUALIZE DAMAGE #
###############################
# Print damage count
cell_damagecount = damage.getDamageCount()
damage.printDamageCount()
# Show complexity distribution
cell_complexitycount = damage.getComplexityDistribution(plot=False)
# Create temp directory if it doesn't exist
temp_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'temp')
os.makedirs(temp_dir, exist_ok=True)

# Produce images of damage. Show to see it instead of saving it. Onlyz in 2D to see only z-projection
#damage.produce3DImage(show=False, saveFile=os.path.join(temp_dir, "xray250_3D.png"), title='X-ray 250 keV - Dose: ' + str(dose) + ' Gy - DSB: ' + str(damage.numDSB))
#damage.produce2DImages(saveFile=os.path.join(temp_dir, "xray250_2D.png"), onlyz=True, title='X-ray 250 keV - Dose: ' + str(dose) + ' Gy - DSB: ' + str(damage.numDSB))
# Get dose-response curves. 'BD' for base damage, 'SSB' for single strand breaks, 'DSB' for double strand breaks
# To see this, recalculatePerEachTrack has to be True (takes longer as the damage is recomputed track by track as the
# dose increases). Skip recal
#damage.getDoseResponseCurve(q='DSB')

# Joining damage in a single SDD file
sddfile = os.path.join(temp_dir, 'test.sdd')
damage.writeSDD(sddfile)
