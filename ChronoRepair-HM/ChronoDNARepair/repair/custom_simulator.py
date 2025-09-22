#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import random
import numpy as np
from ChronoDNARepair.repair.running import Simulator as BaseSimulator
from ChronoDNARepair.repair import output
from ChronoDNARepair.induction.damage import DamageToDNA

class CustomSimulator(BaseSimulator):
    """Custom simulator that handles cellX/runY directory structure.
    
    This simulator version expects a directory structure where damage data is organized as:
    basepath/cell1/run1/{damage files}
    basepath/cell1/run2/{damage files}
    ...
    basepath/cell2/run1/{damage files}
    ...
    
    It allows simulating repair for each cell based on its specific damage files.
    """
    
    def _get_cell_dirs(self, path):
        """Get directories that represent cells (named cellX)"""
        cell_dirs = []
        try:
            for item in os.listdir(path):
                if item.startswith('cell') and os.path.isdir(os.path.join(path, item)):
                    cell_dirs.append(item)
        except Exception as e:
            print(f"Error reading directory {path}: {e}")
        
        # Sort cell directories numerically (cell1, cell2, etc.)
        def extract_cell_num(cell_dir):
            try:
                return int(cell_dir[4:])  # Extract number after 'cell'
            except ValueError:
                return 0
                
        return sorted(cell_dirs, key=extract_cell_num)
    
    def _get_run_dirs(self, cell_path):
        """Get all run directories within a cell directory that contain valid damage data"""
        run_dirs = []
        try:
            for item in os.listdir(cell_path):
                if item.startswith('run'):  # Filter for run directories
                    run_path = os.path.join(cell_path, item)
                    if os.path.isdir(run_path):
                        # Check if directory contains required files
                        damage_file = os.path.join(run_path, 'DNADamage_sdd.txt')
                        phsp_file = os.path.join(run_path, 'DNADamage.phsp')
                        if os.path.exists(damage_file) and os.path.exists(phsp_file):
                            if os.path.getsize(damage_file) > 0:  # only those with actual data
                                run_dirs.append(item)
                else:
                    # If naming convention is different (without 'run' prefix)
                    run_path = os.path.join(cell_path, item)
                    if os.path.isdir(run_path):
                        damage_file = os.path.join(run_path, 'DNADamage_sdd.txt')
                        phsp_file = os.path.join(run_path, 'DNADamage.phsp')
                        if os.path.exists(damage_file) and os.path.exists(phsp_file):
                            if os.path.getsize(damage_file) > 0:
                                run_dirs.append(item)
        except Exception as e:
            print(f"Error reading directory {cell_path}: {e}")
        
        # Sort run directories - try to extract numerical part for proper sorting
        def extract_run_num(run_dir):
            if run_dir.startswith('run'):
                try:
                    return int(run_dir[3:])  # Extract number after 'run'
                except ValueError:
                    pass
            try:
                # Try to extract any numerical part
                return int(''.join(filter(str.isdigit, run_dir)))
            except ValueError:
                return 0
                
        return sorted(run_dirs, key=extract_run_num)
        
    def ReadDamage(self, basepath, maxDose=2.0, version='2.0', recalculatePerEachTrack=False, cell_idx=None):
        """Read damage from the given base path with cell/run structure.
        
        Args:
            basepath: Base directory containing cell directories
            maxDose: Maximum dose to read (can be a single value or a list for multiple radiation types)
            version: SDD format version
            recalculatePerEachTrack: Whether to recalculate damage for each track
            cell_idx: Specific cell index to read (None means use any available cell)
        """
        damage = DamageToDNA(messages=self.messages)
        
        # Handle path and dose as lists or single values
        npaths = 1
        if isinstance(basepath, list):
            npaths = len(basepath)
            if not isinstance(maxDose, list) or len(maxDose) != npaths:
                print('Doses have to match the number of radiations. The same dose will be used for all')
                if not isinstance(maxDose, list):
                    maxDose = [maxDose] * npaths
                else:
                    maxDose = [maxDose[0]] * npaths
        else:
            basepath = [basepath]
            maxDose = [maxDose]

        totalaccumulateddose = 0
        
        for ib, bpath in enumerate(basepath):
            # Get cell directories
            cell_dirs = self._get_cell_dirs(bpath)
            
            if not cell_dirs:
                raise ValueError(f"No cell directories found in {bpath}")
            
            # Use specific cell if requested, otherwise use first available
            if cell_idx is not None:
                target_cell = f"cell{cell_idx}" 
                if target_cell in cell_dirs:
                    cell_dirs = [target_cell]
                else:
                    raise ValueError(f"Cell {target_cell} not found in {bpath}")
            else:
                # Just use the first cell found
                cell_dirs = [cell_dirs[0]]
            
            # Process each cell directory
            for cell_dir in cell_dirs:
                cell_path = os.path.join(bpath, cell_dir)
                run_dirs = self._get_run_dirs(cell_path)
                
                if not run_dirs:
                    print(f"Warning: No valid run directories found in {cell_path}")
                    continue
                    
                print(f"Reading damage from {cell_dir} with {len(run_dirs)} run directories")
                
                # Process each run in the cell
                accumulatedose = 0
                for run_dir in run_dirs:
                    run_path = os.path.join(cell_path, run_dir)
                    
                    # Calculate time for this dose contribution based on accumulated dose so far
                    time = self._getTimeForDose(accumulatedose)
                    if 0 < self.irradiationTime < time:
                        time = 1e20
                    
                    # Read damage data from this run
                    damage.readSDDAndDose(run_path, version=version, particleTime=time, lesionTime=time)
                    
                    # Update accumulated dose
                    accumulatedose = damage.accumulateDose - totalaccumulateddose
                    if maxDose[ib] > 0 and accumulatedose > maxDose[ib]:
                        break
            
            totalaccumulateddose += accumulatedose
        
        # Populate damages
        damage.populateDamages(getVideo=False, stopAtDose=-1, stopAtTime=0, recalculatePerEachTrack=recalculatePerEachTrack)
        self.runManager.damage = damage
        return damage
    
    def Run(self, nCells, rereadDamageForNewRuns=True, basepath=None, maxDose=-1, version='2.0', verbose=1, getVideo=False):
        """
        Run the simulation for multiple cells.
        
        In this custom version, if rereadDamageForNewRuns is True and basepath contains a cell-based structure,
        we will try to use different cell directories for each simulation to better represent cell variations.
        """
        # Check if we can use multiple cells from the cell structure
        if rereadDamageForNewRuns and basepath is not None:
            # Get available cell directories
            if isinstance(basepath, list):
                cell_dirs = self._get_cell_dirs(basepath[0])
            else:
                cell_dirs = self._get_cell_dirs(basepath)
            
            # If we have multiple cell directories available, use them for simulation
            if len(cell_dirs) >= 1:
                # Determine how many cells we can simulate
                available_cells = min(len(cell_dirs), nCells)
                
                if available_cells < nCells:
                    print(f"Warning: Found {available_cells} cell directories, but {nCells} were requested.")
                    print(f"Will simulate {available_cells} cells with different damage data, and reuse some if needed.")
                
                # Initialize necessary properties for output storage
                self.nRuns = nCells  # Store nCells at the simulator level
                self.runManager.TotalRuns = nCells  # Set total runs for progress display
                self.runManager.runoutputDSB = output.AverageTimeCurveOverRuns()
                self.runManager.runoutputFoci = output.AverageTimeCurveOverRuns()
                self.runManager.runoutputMisrepairedDSB = output.AverageTimeCurveOverRuns()
                
                # Initially set plotflag to False - we'll only enable it for the last cell
                self.runManager.plotflag = False
                
                # Track how many cells we've simulated
                simulated_cells = 0
                
                # First simulate with actual cell directories
                for i in range(available_cells):
                    cell_idx = int(cell_dirs[i][4:])  # Extract number from cellX
                    simulated_cells += 1
                    if verbose > 0:
                        print(f"Running cell {simulated_cells}/{nCells} using damage from {cell_dirs[i]}")
                    
                    # Read damage from specific cell directory
                    self.ReadDamage(basepath, maxDose, version, False, cell_idx)
                    
                    # Set runManager to simulate only 1 cell
                    self.runManager.nRuns = 1
                    self.runManager.currentrun = simulated_cells - 1
                    
                    # Only plot for the last cell (when we've reached nCells)
                    plot_this_cell = (simulated_cells == nCells)
                   # self.runManager.plotflag = plot_this_cell
                    
                    # Run for this cell - note the proper capitalization of 'Run'
                    self.runManager.Run(verbose=verbose, 
                                       outputnorm=True,
                                       getVideo=getVideo)
                    
                # If more cells needed than directories available, reread some cells
                for i in range(available_cells, nCells):
                    # Reuse cells in a cycling manner
                    cell_idx = int(cell_dirs[i % available_cells][4:])
                    simulated_cells += 1
                    if verbose > 0:
                        print(f"Running cell {simulated_cells}/{nCells} using reused damage from {cell_dirs[i % available_cells]}")
                    
                    # Read damage from specific cell directory
                    self.ReadDamage(basepath, maxDose, version, False, cell_idx)
                    
                    # Set runManager to simulate only 1 cell
                    self.runManager.nRuns = 1
                    self.runManager.currentrun = simulated_cells - 1
                    
                    # Only plot for the last cell (when we've reached nCells)
                    plot_this_cell = (simulated_cells == nCells)
                    self.runManager.plotflag = plot_this_cell
                    
                    # Run for this cell - note the proper capitalization of 'Run'
                    self.runManager.Run(verbose=verbose,
                                       outputnorm=True,
                                       getVideo=getVideo)
                
                # Set the final output attributes so they can be accessed from the main script
                self.avgRemainingDSBOverTime = self.runManager.runoutputDSB
                self.avgRemainFociOverTime = self.runManager.runoutputFoci
                self.avgMisrepairedDSBOverTime = self.runManager.runoutputMisrepairedDSB
                self.celloutput = self.runManager.cellcultureoutput
                
                return
        
        # If we can't use the cell structure, fall back to the original behavior
        super().Run(nCells, rereadDamageForNewRuns, basepath, maxDose, version, verbose, getVideo)