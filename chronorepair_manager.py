#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChronoRepair Manager Module

This module provides functions to interact with ChronoDNARepair framework for DNA repair simulations.
It is designed to be used with TOPAS simulation data organized in cell/run directories.

Functions:
- setup_repair_simulation: Set up parameters for DNA repair simulation
- run_repair_simulation: Run DNA repair simulations and return results
- process_cell_repair: Process repair simulation for a single cell
- process_multicell_repair: Process repair simulations for multiple cells
- compute_survival_metrics: Compute cell survival metrics from repair results
- create_repair_visualizations: Create visualizations of repair process

Examples:
    # Set up repair simulation parameters
    sim_params = setup_repair_simulation(
        exposure_time=23,
        nucleus_max_radius=4.65,
        initial_dose_rate=0.13803/3600,
        half_life=(59.39*24)*3600
    )
    
    # Run repair simulation
    repair_results = run_repair_simulation(
        sim_params=sim_params,
        damage_path='/path/to/damage/files',
        n_cells=10
    )
    
    # Display results
    display_repair_results(repair_results)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict

# Add ChronoRepair-HM paths
current_dir = os.path.dirname(os.path.abspath(__file__))
chrono_dir = os.path.join(current_dir, 'ChronoRepair-HM')
sys.path.insert(0, chrono_dir)
chrono_dna_dir = os.path.join(chrono_dir, 'ChronoDNARepair')
sys.path.insert(0, chrono_dna_dir)

try:
    # Import necessary ChronoDNARepair modules
    from ChronoDNARepair.repair.running import Simulator as BaseSimulator
    from ChronoDNARepair.repair import output
    from ChronoDNARepair.induction.damage import DamageToDNA
    
    # Import custom simulator
    sys.path.insert(0, os.path.join(chrono_dir, 'examples'))
    from custom_simulator import CustomSimulator
except ImportError as e:
    print(f"Warning: ChronoDNARepair module import error: {e}")
    print("Make sure ChronoRepair-HM is properly installed.")


def setup_repair_simulation(
    exposure_time=23,              # Hours
    simulation_time=23,            # Hours (usually matches exposure_time)
    time_steps=23,                 # Number of time steps for simulation
    nucleus_max_radius=4.65,       # Microns
    diffusion_model='free',        # Model for DNA fragment diffusion
    dsb_model='standard',          # Model for DSB repair
    ssb_model='standard',          # Model for SSB repair
    bd_model='standard',           # Model for base damage repair
    dose_rate_function='exponential',  # Function to model dose rate over time
    initial_dose_rate=0.13803/3600,    # Gy/hour (here converted to Gy/second)
    half_life=(59.39*24)*3600      # Half-life in seconds (here 59.39 days)
):
    """
    Set up parameters for a DNA repair simulation.
    
    Args:
        exposure_time: Duration of radiation exposure in hours
        simulation_time: Total simulation time in hours (often equal to exposure_time)
        time_steps: Number of discrete time steps for the simulation
        nucleus_max_radius: Maximum radius of cell nucleus in microns
        diffusion_model: Model for DNA fragment diffusion ('free', 'confined', etc.)
        dsb_model: Model for double-strand break repair ('standard', 'fast', etc.)
        ssb_model: Model for single-strand break repair ('standard', 'fast', etc.)
        bd_model: Model for base damage repair ('standard', 'fast', etc.)
        dose_rate_function: Function to model dose rate over time ('constant', 'exponential', etc.)
        initial_dose_rate: Initial dose rate in Gy/hour
        half_life: Half-life of radioactive source in seconds
        
    Returns:
        Dictionary containing all simulation parameters
    """
    # Convert times to seconds
    initial_time = 0
    final_time = simulation_time * 3600  # Convert hours to seconds
    irradiation_time = exposure_time * 3600  # Convert hours to seconds
    
    # Combine parameters in a dictionary
    params = {
        'time_options': [initial_time, final_time, time_steps],
        'nucleus_max_radius': nucleus_max_radius,
        'diffusion_model': diffusion_model,
        'dsb_model': dsb_model,
        'ssb_model': ssb_model,
        'bd_model': bd_model,
        'dose_rate_function': dose_rate_function,
        'initial_dose_rate': initial_dose_rate,  # Gy/s
        'half_life': half_life,  # seconds
        'irradiation_time': irradiation_time  # seconds
    }
    
    return params


def run_repair_simulation(
    sim_params,            # Parameters from setup_repair_simulation
    damage_path,           # Path to directory with cell/run structure
    n_cells=10,            # Number of cells to simulate
    maximum_dose=-1,       # Maximum dose to consider (-1 for all)
    verbose=1,             # Verbosity level
    get_video=False        # Generate visualization video
):
    """
    Run DNA repair simulation for multiple cells.
    
    Args:
        sim_params: Dictionary of simulation parameters from setup_repair_simulation
        damage_path: Path to damage files with cell/run structure
        n_cells: Number of cells to simulate
        maximum_dose: Maximum dose to consider (-1 for all)
        verbose: Verbosity level (0=quiet, 1=normal, 2=detailed)
        get_video: Generate visualization video
    
    Returns:
        Dictionary containing simulation results
    """
    # Initialize simulator with parameters
    sim = CustomSimulator(
        timeOptions=sim_params['time_options'],
        diffusionmodel=sim_params['diffusion_model'],
        dsbmodel=sim_params['dsb_model'],
        ssbmodel=sim_params['ssb_model'],
        bdmodel=sim_params['bd_model'],
        nucleusMaxRadius=sim_params['nucleus_max_radius'],
        doseratefunction=sim_params['dose_rate_function'],
        doseratefunctionargs=[sim_params['initial_dose_rate'], sim_params['half_life']],
        irradiationTime=sim_params['irradiation_time']
    )
    
    # Run the simulation
    sim.Run(
        nCells=n_cells,
        rereadDamageForNewRuns=True,  # Use different damage data for each cell when available
        basepath=damage_path,
        maxDose=maximum_dose,
        verbose=verbose,
        getVideo=get_video
    )
    
    # Collect results
    results = {
        'dsb_output': sim.avgRemainingDSBOverTime,
        'foci_output': sim.avgRemainFociOverTime,
        'misrepaired_output': sim.avgMisrepairedDSBOverTime,
        'cell_output': sim.celloutput,
        'parameters': sim_params
    }
    
    return results


def process_cell_repair(damage_path, sim_params, maximum_dose=-1, verbose=0):
    """
    Process repair simulation for a single cell.
    
    Args:
        damage_path: Path to cell directory containing run directories with damage files
        sim_params: Dictionary of simulation parameters
        maximum_dose: Maximum dose to consider (-1 for all)
        verbose: Verbosity level
    
    Returns:
        Dictionary with repair results for the cell
    """
    # Run repair simulation for a single cell
    results = run_repair_simulation(
        sim_params=sim_params,
        damage_path=damage_path,
        n_cells=1,
        maximum_dose=maximum_dose,
        verbose=verbose,
        get_video=False
    )
    
    # Extract DSB remaining data
    dsb_output = results['dsb_output']
    times = dsb_output.times
    avg_dsb_remaining = dsb_output.avgyvalues
    
    # Compute cell survival
    cell_output = results['cell_output']
    surviving = cell_output.celllist[0].Surviving if len(cell_output.celllist) > 0 else False
    
    # Combine results into simple dictionary
    cell_results = {
        'times': times,
        'dsb_remaining': avg_dsb_remaining,
        'surviving': surviving,
        'foci': results['foci_output'].avgyvalues if hasattr(results['foci_output'], 'avgyvalues') else None,
        'misrepaired': results['misrepaired_output'].avgyvalues if hasattr(results['misrepaired_output'], 'avgyvalues') else None,
    }
    
    return cell_results


def process_multicell_repair(basepath, n_cells, sim_params, maximum_dose=-1, verbose=1):
    """
    Process repair simulations for multiple cells.
    
    Args:
        basepath: Path to directory containing cell directories
        n_cells: Number of cells to process
        sim_params: Dictionary of simulation parameters
        maximum_dose: Maximum dose to consider (-1 for all)
        verbose: Verbosity level
        
    Returns:
        Dictionary with aggregated repair results for all cells
    """
    # Run repair simulation for multiple cells
    results = run_repair_simulation(
        sim_params=sim_params,
        damage_path=basepath,
        n_cells=n_cells,
        maximum_dose=maximum_dose,
        verbose=verbose,
        get_video=False
    )
    
    # Extract DSB remaining data
    dsb_output = results['dsb_output']
    times = dsb_output.times
    avg_dsb_remaining = dsb_output.avgyvalues
    
    # Calculate survival statistics
    cell_output = results['cell_output']
    surviving_cells = sum(1 for cell in cell_output.celllist if cell.Surviving)
    total_cells = len(cell_output.celllist)
    survival_fraction = surviving_cells / total_cells if total_cells > 0 else 0
    
    # Compile results
    multicell_results = {
        'times': times,
        'dsb_remaining': avg_dsb_remaining,
        'foci': results['foci_output'].avgyvalues if hasattr(results['foci_output'], 'avgyvalues') else None,
        'misrepaired': results['misrepaired_output'].avgyvalues if hasattr(results['misrepaired_output'], 'avgyvalues') else None,
        'survival': {
            'surviving_cells': surviving_cells,
            'total_cells': total_cells,
            'survival_fraction': survival_fraction
        }
    }
    
    return multicell_results


def compute_survival_metrics(repair_results):
    """
    Compute cell survival metrics from repair results.
    
    Args:
        repair_results: Results from repair simulation
        
    Returns:
        Dictionary with survival metrics
    """
    cell_output = repair_results['cell_output']
    surviving_cells = sum(1 for cell in cell_output.celllist if cell.Surviving)
    total_cells = len(cell_output.celllist)
    survival_fraction = surviving_cells / total_cells if total_cells > 0 else 0
    
    survival_metrics = {
        'surviving_cells': surviving_cells,
        'total_cells': total_cells,
        'survival_fraction': survival_fraction,
        'survival_percent': survival_fraction * 100
    }
    
    return survival_metrics


def display_repair_results(repair_results, title=None):
    """
    Display comprehensive results from repair simulations.
    
    Args:
        repair_results: Results from repair simulation
        title: Optional title for the output
    """
    if title:
        print(f"\n{title}")
        print("="* len(title))
    
    # DSB repair data
    if 'dsb_output' in repair_results:
        dsb_output = repair_results['dsb_output']
        times = dsb_output.times
        avg_dsb_remaining = dsb_output.avgyvalues
        
        print("\nDSB Repair Over Time:")
        print("-" * 40)
        for t in range(len(times)):
            # Convert time to hours for display
            time_h = times[t] / 3600 if times[t] > 100 else times[t]  # If already in hours, don't convert
            time_unit = 'h' if times[t] > 100 else 'h'  # Assume times are in seconds if large
            print(f'Time: {time_h:.2f} {time_unit}, Fraction of DSB remaining: {avg_dsb_remaining[t]:.4f}')
    
    # Calculate survival statistics
    survival_metrics = compute_survival_metrics(repair_results)
    
    print("\n----- Cell Survival Results -----")
    print(f"Surviving cells: {survival_metrics['surviving_cells']}/{survival_metrics['total_cells']}")
    print(f"Survival fraction: {survival_metrics['survival_fraction']:.4f} ({survival_metrics['survival_percent']:.1f}%)")
    print("--------------------------------")


def plot_dsb_repair_kinetics(repair_results, title='DSB Repair Kinetics'):
    """
    Plot DSB repair kinetics from repair results.
    
    Args:
        repair_results: Results from repair simulation
        title: Title for the plot
    """
    if 'dsb_output' not in repair_results:
        print("No DSB repair data available to plot.")
        return
    
    dsb_output = repair_results['dsb_output']
    times = dsb_output.times
    avg_dsb_remaining = dsb_output.avgyvalues
    
    # Convert times to hours if needed
    times_h = [t/3600 if t > 100 else t for t in times]  # Convert to hours if in seconds
    
    plt.figure(figsize=(10, 6))
    plt.plot(times_h, avg_dsb_remaining, 'o-', linewidth=2, markersize=8)
    plt.xlabel('Time (hours)', fontsize=12)
    plt.ylabel('Fraction of DSBs remaining', fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1.05)
    
    # Add text with survival rate
    if 'cell_output' in repair_results:
        survival_metrics = compute_survival_metrics(repair_results)
        sf = survival_metrics['survival_fraction']
        plt.text(0.05, 0.05, f'Survival Fraction: {sf:.4f} ({sf*100:.1f}%)', 
                transform=plt.gca().transAxes, fontsize=12,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray'))
    
    plt.tight_layout()
    plt.show()


def compare_repair_scenarios(scenarios_results, labels=None, title='Comparison of DSB Repair Kinetics'):
    """
    Compare DSB repair kinetics between different scenarios.
    
    Args:
        scenarios_results: List of repair results from different scenarios
        labels: List of labels for each scenario
        title: Title for the plot
    """
    plt.figure(figsize=(12, 7))
    
    # Prepare labels if not provided
    if labels is None:
        labels = [f'Scenario {i+1}' for i in range(len(scenarios_results))]
    
    # Prepare colors
    colors = plt.cm.viridis(np.linspace(0, 1, len(scenarios_results)))
    
    # Plot each scenario
    for i, results in enumerate(scenarios_results):
        if 'dsb_output' not in results:
            print(f"No DSB repair data available for {labels[i]}.")
            continue
        
        dsb_output = results['dsb_output']
        times = dsb_output.times
        avg_dsb_remaining = dsb_output.avgyvalues
        
        # Convert times to hours if needed
        times_h = [t/3600 if t > 100 else t for t in times]  # Convert to hours if in seconds
        
        # Plot with legend
        plt.plot(times_h, avg_dsb_remaining, 'o-', linewidth=2, markersize=8, 
                 color=colors[i], label=f'{labels[i]} (SF: {compute_survival_metrics(results)["survival_fraction"]:.2f})')
    
    plt.xlabel('Time (hours)', fontsize=12)
    plt.ylabel('Fraction of DSBs remaining', fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1.05)
    plt.legend(fontsize=10, loc='upper right')
    
    plt.tight_layout()
    plt.show()


def calculate_repair_half_time(repair_results):
    """
    Calculate the repair half-time (time to repair 50% of DSBs).
    
    Args:
        repair_results: Results from repair simulation
        
    Returns:
        Half-time in hours
    """
    if 'dsb_output' not in repair_results:
        print("No DSB repair data available to calculate half-time.")
        return None
    
    dsb_output = repair_results['dsb_output']
    times = dsb_output.times
    avg_dsb_remaining = dsb_output.avgyvalues
    
    # Find when remaining DSBs drop below 50%
    half_time = None
    for i in range(len(avg_dsb_remaining)):
        if avg_dsb_remaining[i] <= 0.5:
            if i > 0:
                # Linear interpolation to get a more precise half-time
                t1, t2 = times[i-1], times[i]
                r1, r2 = avg_dsb_remaining[i-1], avg_dsb_remaining[i]
                
                if r1 != r2:  # Avoid division by zero
                    half_time = t1 + (0.5 - r1) * (t2 - t1) / (r2 - r1)
                else:
                    half_time = t1
            else:
                half_time = times[i]
            break
    
    # Convert to hours if in seconds
    if half_time is not None and half_time > 100:
        half_time /= 3600
        
    return half_time


if __name__ == "__main__":
    # Example usage
    print("This module provides functions for ChronoDNARepair simulations.")
    print("To use it, import the module and call its functions.")
