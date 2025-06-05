import os
import numpy as np
import pandas as pd
from chemistry_output_manager import read_irtgvalue_phase_space, filter_inactive_species, create_pivots, get_final_values
from phsp_manager import count_phsp_particles
from topas_csv_files_manager import process_csv_file, process_original_hists
from dnadamage_phsp_manager import read_dnadamage_phase_space, merge_dnadamage_files

def multirun_processing(maxruns, filebase):
    """Process multiple TOPAS simulation runs and aggregate results.
    
    Args:
        maxruns (int): Number of runs to process
        filebase (str): Base directory containing the run directories
    
    Returns:
        dict: Dictionary containing aggregated results from all runs
    """
    # Generate run directories
    runs_filebases = [os.path.join(filebase, f'run{i+1}/') for i in range(maxruns)]
    
    # Species of interest for G-Values
    species_of_interest = ['OH^-1', 'H2O2^0', 'e_aq^-1', 'H^0', 'H3O^1','OH^0', 'H_2^0']
    
    Cell_results_files = {
        'Original_hists': 'DoseToCell_I125Beam.csv',
        'DoseToNucl_ph2': 'DoseNucleus_Total.csv',
        'DoseToNucl_ph3': 'DoseNucleus_Ph3.csv',
        'Ecell': 'EnergyToCell.csv',
        'NP_el': 'PhaseSpace_NP',
        'GValues': 'IRTGValue',
        'DNADamage': 'DNADamage'
    }
    
    # Initialize results structure
    Cell_results = {
        'Original_hists': {'value': 0},  # Single value, no error needed
        'DoseToNucl_ph2': {'value': 0.0, 'error': 0.0},
        'DoseToNucl_ph3': {'value': 0.0, 'error': 0.0},
        'Ecell': {'value': 0.0, 'error': 0.0},
        'NP_el': {'value': 0},  # Count value, no error needed
        'GValues': {species: {'value': 0.0, 'error': 0.0} for species in species_of_interest},
        'DNADamage': {
            'totals': {},  # Will store the total damage statistics
            'frames': [],  # Will store DataFrames from each run
            'damage_cols': [
                'DSBs', 'DSBs_Direct', 'DSBs_Indirect', 'DSBs_Hybrid',
                'SSBs', 'SSBs_Direct', 'SSBs_QuasiDirect', 'SSBs_Indirect',
                'SBs', 'SBs_Direct', 'SBs_QuasiDirect', 'SBs_Indirect',
                'SSB+s', 'DSB+s', 'More complex damages', 
                'BDs', 'BDs_Direct', 'BDs_QuasiDirect', 'BDs_Indirect',
                'Foci_150nm', 'Foci_500nm'
            ]
        }
    }
    
    # Process all runs
    for run_idx, run_base in enumerate(runs_filebases):
        print(f"Processing run {run_idx+1}/{maxruns}...", end='\r')
        
        # Process each result type
        for result_key, file_name in Cell_results_files.items():
            file_path = os.path.join(run_base, file_name)
            
            if result_key == 'GValues':
                try:
                    # Process chemical phase (G-Values)
                    df = read_irtgvalue_phase_space(file_path)
                    df = filter_inactive_species(df)
                    pivot, pivot_err = create_pivots(df)
                    summary, summary_err, last_time = get_final_values(pivot, pivot_err)
                    
                    # Accumulate results for each species
                    for species in species_of_interest:
                        if species in summary.index:
                            if 'values' not in Cell_results['GValues'][species]:
                                Cell_results['GValues'][species]['values'] = []
                                Cell_results['GValues'][species]['errors'] = []
                            Cell_results['GValues'][species]['values'].append(summary.loc[species])
                            Cell_results['GValues'][species]['errors'].append(summary_err.loc[species])
                except Exception as e:
                    print(f"\nError processing GValues for run {run_idx+1}: {e}")
                    continue
            
            elif result_key == 'NP_el':
                # Process phase space files
                count, _ = count_phsp_particles(file_path)
                Cell_results[result_key]['value'] += count
            
            elif result_key == 'DNADamage':
                try:
                    # Process DNA damage files
                    df = read_dnadamage_phase_space(file_path)
                    Cell_results['DNADamage']['frames'].append(df)
                except Exception as e:
                    print(f"\nError processing DNADamage for run {run_idx+1}: {e}")
                    continue
                
            elif result_key == 'Original_hists':
                # Process Original_hists
                Cell_results[result_key]['value'] += process_original_hists(file_path)
            
            else:
                # Process other CSV files (dose and energy)
                stats = process_csv_file(file_path, return_stats_object=True)
                if not stats:
                    continue
                
                # Accumulate values for statistical processing
                if 'values' not in Cell_results[result_key]:
                    Cell_results[result_key]['values'] = []
                    Cell_results[result_key]['stats'] = stats
                else:
                    Cell_results[result_key]['values'].append(stats.mean)
                    Cell_results[result_key]['stats'] += stats
    
    print("\nProcessing complete!")
    
    # Calculate final statistics
    for key in ['DoseToNucl_ph2', 'DoseToNucl_ph3', 'Ecell']:
        if 'values' in Cell_results[key]:
            values = np.array(Cell_results[key]['values'])
            stats = Cell_results[key]['stats']
            
            # Calculate mean and error
            Cell_results[key]['value'] = stats.mean
            std_error = stats.standard_deviation / np.sqrt(stats.histories_with_scorer_active)
            run_variation = np.std(values) if len(values) > 1 else 0
            Cell_results[key]['error'] = 2 * np.sqrt(std_error**2 + run_variation**2)  # 2σ confidence
            
            # Clean up temporary storage
            del Cell_results[key]['values']
            del Cell_results[key]['stats']
    
    # Calculate final G-Values statistics
    for species in species_of_interest:
        if 'values' in Cell_results['GValues'][species]:
            values = np.array(Cell_results['GValues'][species]['values'])
            errors = np.array(Cell_results['GValues'][species]['errors'])
            
            if len(values) > 0:
                # Calculate final statistics
                mean_value = np.mean(values)
                sem = np.std(values) / np.sqrt(len(values))
                weighted_error = np.sqrt(np.sum(errors**2)) / len(errors)
                total_error = np.sqrt(sem**2 + weighted_error**2)
                
                # Store final values
                Cell_results['GValues'][species]['value'] = mean_value
                Cell_results['GValues'][species]['error'] = 2 * total_error  # 2σ confidence
                
                # Clean up temporary storage
                del Cell_results['GValues'][species]['values']
                del Cell_results['GValues'][species]['errors']
    
    # Process DNA damage results if we have collected any frames
    if Cell_results['DNADamage']['frames']:
        # Merge all DNA damage frames
        merged_df = pd.concat(Cell_results['DNADamage']['frames'], ignore_index=True)
        
        # Calculate totals for each damage column
        for col in Cell_results['DNADamage']['damage_cols']:
            if col in merged_df.columns:
                total = merged_df[col].sum()
                Cell_results['DNADamage']['totals'][col] = float(total)
        
        # Calculate energy statistics
        if 'Energy_imparted_per_event [keV]' in merged_df.columns:
            total_energy = merged_df['Energy_imparted_per_event [keV]'].sum() / 1000  # Convert to MeV
            Cell_results['DNADamage']['totals']['energy'] = {
                'total_MeV': float(total_energy)
            }
        
        # Store total number of events
        Cell_results['DNADamage']['totals']['total_events'] = len(merged_df)
        
        # Clean up the frames to save memory
        del Cell_results['DNADamage']['frames']
    
    return Cell_results

def multicell_processing(maxcells, maxruns, filebase):
    """Process multiple cells each containing multiple TOPAS simulation runs and aggregate results.
    
    Args:
        maxcells (int): Number of cell directories to process
        maxruns (int): Number of runs to process per cell
        filebase (str): Base directory containing the cell directories
        
    Returns:
        list: List of dictionaries containing aggregated results for each cell
    """
    # Generate cell directories
    cell_filebases = [os.path.join(filebase, f'cell{i+1}/') for i in range(maxcells)]
    
    # List to store results for each cell
    all_cell_results = []
    
    # Process runs for each cell
    for cell_idx, cell_base in enumerate(cell_filebases):
        print(f"\nProcessing cell {cell_idx+1}/{maxcells}...")
        try:
            # Call multirun_processing for each cell
            cell_results = multirun_processing(maxruns, cell_base)
            all_cell_results.append(cell_results)
        except Exception as e:
            print(f"Error processing cell {cell_idx+1}: {e}")
            continue
            
    print("\nAll cells processed!")
    return all_cell_results

def process_multicell_results(all_cell_results):
    """Process and aggregate results from multiple cells to compute statistics.
    
    Args:
        all_cell_results (list): List of Cell_results dictionaries from multicell_processing
        
    Returns:
        dict: Dictionary containing mean and standard deviation for each quantity across cells
    """
    # Initialize aggregated results structure
    aggregated_results = {
        'Original_hists': {'mean': 0.0, 'std_dev': 0.0},
        'DoseToNucl_ph2': {'mean': 0.0, 'std_dev': 0.0},
        'DoseToNucl_ph3': {'mean': 0.0, 'std_dev': 0.0},
        'Ecell': {'mean': 0.0, 'std_dev': 0.0},
        'NP_el': {'mean': 0.0, 'std_dev': 0.0},
        'GValues': {},
        'DNADamage': {
            'totals': {}
        }
    }
    
    # Get number of cells
    n_cells = len(all_cell_results)
    if n_cells == 0:
        return aggregated_results
        
    # First, collect all values for each quantity
    collected_values = {
        'Original_hists': [],
        'DoseToNucl_ph2': [],
        'DoseToNucl_ph3': [],
        'Ecell': [],
        'NP_el': [],
        'GValues': {},
        'DNADamage': {
            'totals': {}
        }
    }
    
    # Initialize GValues species list from first cell
    species_list = list(all_cell_results[0]['GValues'].keys())
    for species in species_list:
        collected_values['GValues'][species] = {'values': [], 'errors': []}
        aggregated_results['GValues'][species] = {'mean': 0.0, 'std_dev': 0.0, 'error': 0.0}
    
    # Initialize DNA damage columns from first cell
    damage_cols = all_cell_results[0]['DNADamage']['totals'].keys()
    for col in damage_cols:
        collected_values['DNADamage']['totals'][col] = []
        aggregated_results['DNADamage']['totals'][col] = {'mean': 0.0, 'std_dev': 0.0}
    
    # Collect values from each cell
    for cell_results in all_cell_results:
        # Simple quantities
        for key in ['Original_hists', 'DoseToNucl_ph2', 'DoseToNucl_ph3', 'Ecell', 'NP_el']:
            collected_values[key].append(cell_results[key]['value'])
        
        # G-Values
        for species in species_list:
            collected_values['GValues'][species]['values'].append(cell_results['GValues'][species]['value'])
            collected_values['GValues'][species]['errors'].append(cell_results['GValues'][species]['error'])
        
        # DNA Damage
        for col in damage_cols:
            if col in cell_results['DNADamage']['totals']:
                value = cell_results['DNADamage']['totals'][col]
                # Handle nested energy dictionary
                if isinstance(value, dict):
                    value = value['total_MeV']
                collected_values['DNADamage']['totals'][col].append(value)
    
    # Calculate statistics
    # Simple quantities
    for key in ['Original_hists', 'DoseToNucl_ph2', 'DoseToNucl_ph3', 'Ecell', 'NP_el']:
        values = np.array(collected_values[key])
        aggregated_results[key]['mean'] = np.mean(values)
        aggregated_results[key]['std_dev'] = np.std(values, ddof=1) if len(values) > 1 else 0
    
    # G-Values
    for species in species_list:
        values = np.array(collected_values['GValues'][species]['values'])
        errors = np.array(collected_values['GValues'][species]['errors'])
        
        # Calculate mean and standard deviation across cells
        mean_value = np.mean(values)
        std_dev = np.std(values, ddof=1) if len(values) > 1 else 0
        
        # Propagate errors
        total_error = np.sqrt(np.mean(errors**2) + std_dev**2)
        
        aggregated_results['GValues'][species]['mean'] = mean_value
        aggregated_results['GValues'][species]['std_dev'] = std_dev
        aggregated_results['GValues'][species]['error'] = total_error
    
    # DNA Damage
    for col in damage_cols:
        values = np.array(collected_values['DNADamage']['totals'][col])
        aggregated_results['DNADamage']['totals'][col]['mean'] = np.mean(values)
        aggregated_results['DNADamage']['totals'][col]['std_dev'] = np.std(values, ddof=1) if len(values) > 1 else 0
    
    return aggregated_results
