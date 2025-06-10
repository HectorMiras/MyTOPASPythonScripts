import pandas as pd
import numpy as np
import re
import os
import random
import sys
import matplotlib.pyplot as plt


# Add ChronoRepair-HM paths
current_dir = os.path.dirname(os.path.abspath(__file__))
chrono_dir = os.path.join(current_dir, 'ChronoRepair-HM')
sys.path.insert(0, chrono_dir)
chrono_dna_dir = os.path.join(chrono_dir, 'ChronoDNARepair')
sys.path.insert(0, chrono_dna_dir)

try:
    from ChronoDNARepair.induction.damage import DamageToDNA
except ImportError:
    print("Warning: ChronoDNARepair module not found. Make sure ChronoRepair-HM is properly installed.")



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

def process_multirun_dnadamage(damagepath='', maxruns=0, dose=-1, random_sampling=True):
    """
    Process multi-run DNA damage simulations from SDD files
    
    Args:
        damagepath: Path to the parent directory containing run subdirectories
        maxruns: Maximum number of runs to process (default=100, use -1 for all available)
        dose: Target dose to stop at (-1 for all available data)
        random_sampling: Whether to sample runs randomly (default=True)
        
    Returns:
        Dictionary containing aggregated damage statistics
    """
    # Define colormap for visualization
    MyColorMap = 'viridis'

    # Initialize damage object
    damage = DamageToDNA()
    
    # Read SDD files - get directories that contain both dose and SDD files
    listOfAvailableDirs = []
    for dir_name in get_numbered_dirs(damagepath):
        newpath = os.path.join(damagepath, dir_name)
        files = os.listdir(newpath)
        if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
            if os.path.getsize(os.path.join(newpath, 'DNADamage_sdd.txt')) > 0:  # only those with actual data
                listOfAvailableDirs.append(dir_name)
    
    # Process up to maxruns directories
    if len(listOfAvailableDirs) == 0:
        print("No valid directories found.")
        return None
    
    # Determine which directories to process
    if random_sampling:
        # Randomize order of directories
        neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
    else:
        # Keep original order
        neworder = listOfAvailableDirs.copy()
    
    # Limit number of runs based on maxruns parameter
    if maxruns > 0 and maxruns < len(neworder):
        neworder = neworder[:maxruns]
        print(f"Processing {maxruns} directories (limited by maxruns).")
    else:
        print(f"Processing all {len(neworder)} available directories.")
        
    # Process each directory
    print(f"Number of directories with damage data to process: {len(neworder)}")
    for i, e in enumerate(neworder):
        path = os.path.join(damagepath, e) + '/'
        #print(f"Reading run {i+1}/{len(neworder)}: {e}")
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
    #damage.printDamageCount()
    # Show complexity distribution
    cell_complexitydist= damage.getComplexityDistribution(plot=False)
    for key, value in cell_complexitydist.items():
        cell_damagecount[f'Complexity{key}'] = value


    # Joining damage in a single SDD file
    sddfile = os.path.join(damagepath, 'test.sdd')
    damage.writeSDD(sddfile)

    return cell_damagecount