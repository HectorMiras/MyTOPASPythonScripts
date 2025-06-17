"""
This module contains display functions for plotting and formatting cell simulation results.
Functions are used to create plots and tables for showing:
- Basic simulation results
- DNA damage distribution
- G-values
- Chemical species across cells
- Enhancement ratios and comparisons
"""

import sys
import os
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.lines import Line2D
#from IPython.display import display
import re
from collections import defaultdict

# Set default figure size and style for better display
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['figure.dpi'] = 100
plt.style.use('seaborn-v0_8-whitegrid')  # Modern style for better visualization

# Avoid warnings about too many open figures
plt.rcParams['figure.max_open_warning'] = 0

def get_output_dirs(output_path=None):
    """Get the paths for plots and tables directories.
    
    Args:
        output_path: Optional path where to create plots and tables directories.
                    If None, creates them in the current directory.
    
    Returns:
        Tuple of (plots_dir, tables_dir)
    """
    if output_path:
        base_dir = output_path
    else:
        base_dir = os.getcwd()
    
    plots_dir = os.path.join(base_dir, 'plots')
    tables_dir = os.path.join(base_dir, 'tables')
    
    # Create directories if they don't exist
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    
    return plots_dir, tables_dir

# Determine if we're running in Jupyter or terminal
def is_jupyter():
    """Check if the code is running in a Jupyter notebook."""
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            return True
        return False
    except ImportError:
        return False

def save_or_show_plot(fig, plot_name, output_path=None):
    """Save the plot to a file and display it in terminal, or just display it in Jupyter.
    
    Args:
        fig: The matplotlib figure to save/show
        plot_name: The name for the plot file
        output_path: Optional path where to save the plot. If None, uses current directory.
    """
    if not is_jupyter():
        plots_dir, _ = get_output_dirs(output_path)
        filename = os.path.join(plots_dir, f'{plot_name}.png')
        fig.savefig(filename)
        print(f"Plot saved to {filename}")
        plt.show()  # Show the plot in terminal
        plt.close(fig)  # Close the figure to free memory
    else:
        plt.show()

def display_results(results):
    """Display processed results with proper formatting and error reporting.
    
    Args:
        results: Dictionary containing simulation results with the following structure:
            - Original_hists: Dict with 'value' key
            - NP_el: Dict with 'value' key
            - DoseToNucl_ph2, DoseToNucl_ph3, Ecell: Dict with 'value' and 'error' keys
            - GValues: Dict of chemical species with 'value' and 'error' keys
            - DNADamage: Optional dict with damage statistics
    """
    print("\nResults Summary:")
    print("-" * 50)
    
    # Display physical quantities
    print("\nPhysical Quantities:")
    print("-" * 50)
    
    # Display particle counts
    print(f"Original histories: {results['Original_hists']['value']:,}")
    print(f"Nanoparticle electrons: {results['NP_el']['value']:,} particles")
    
    # Display dose and energy measurements
    for key in ['DoseToNucl_ph2', 'DoseToNucl_ph3', 'Ecell']:
        print(f"\n{key}:")
        print(f"  - Value: {results[key]['value']:.6e} ± {results[key]['error']:.6e} (2σ)")
    
    # Display chemical phase results 
    print("\nChemical Phase Results (G-Values):")
    print("-" * 50)
    
    # Sort species by G-Value for better presentation
    species_data = [(species, data['value'], data['error']) 
                   for species, data in results['GValues'].items() 
                   if 'value' in data]
    species_data.sort(key=lambda x: x[1], reverse=True)
    
    for species, value, error in species_data:
        print(f"\n{species}:")
        print(f"  - G-Value: {value:.4f} ± {error:.4f} molecules/100eV (2σ)")
    
    # Display DNA damage results if available
    if 'DNADamage' in results:
        print("\nDNA Damage Results:")
        print("-" * 50)
        
        if 'Dose' in results['DNADamage']:
            dnadose = results['DNADamage']['Dose']
            print(f"Total dose deposited: {dnadose:.2f} Gy")
        
        # Group and display damage statistics
        damage_keys = ['DSB', 'DSB_Direct', 'DSB_Indirect', 'DSB_Hybrid', 'SSB', 'SSB_Direct', 'SSB_Indirect', 
                      'SB', 'SB_Direct', 'SB_Indirect', 'BD', 'BD_Direct', 'BD_Indirect', 'DSB_positions', 
                      'Number_of_foci', 'Complexity2', 'Complexity3', 'Complexity4', 'Complexity5', 'Complexity6',
                      'Complexity7', 'Complexity8', 'Complexity9', 'Complexity10', 'Complexity11', 'Complexity12', 
                      'Complexity13', 'Complexity14', 'Complexity15']
        
        for dmg in damage_keys:
            if dmg in results['DNADamage'].keys():
                val = results['DNADamage'][dmg]
                print(f"\n{dmg}: {val:.0f}")

def plot_damage_distribution(damage_totals, save_plots=False):
    """Create a stacked bar plot showing direct vs indirect damage distribution.
    
    Args:
        damage_totals: Dictionary containing damage statistics with keys like 
                      'DSB_Direct', 'DSB_Indirect', etc.
        save_plots: If True, saves plots to files instead of displaying them.
    """
    damage_pairs = [
        ('DSB_Direct', 'DSB_Indirect'),
        ('SSB_Direct', 'SSB_Indirect'),
        ('SB_Direct', 'SB_Indirect'),
        ('BD_Direct', 'BD_Indirect')
    ]
    
    valid_pairs = [(direct, indirect) for direct, indirect in damage_pairs 
                  if direct in damage_totals and indirect in damage_totals]
    
    if valid_pairs:
        fig, ax = plt.subplots(figsize=(10, 6))
        bar_width = 0.35
        x = np.arange(len(valid_pairs))
        labels = [pair[0].split('_')[0] for pair in valid_pairs]
        
        # Plot stacked bars with consistent colors
        direct_color = '#1f77b4'  # blue
        indirect_color = '#ff7f0e'  # orange
        for i, (direct, indirect) in enumerate(valid_pairs):
            direct_sum = damage_totals[direct]
            indirect_sum = damage_totals[indirect]
            ax.bar(i, direct_sum, bar_width, color=direct_color, label='Direct' if i == 0 else "")
            ax.bar(i, indirect_sum, bar_width, bottom=direct_sum, color=indirect_color, 
                  label='Indirect' if i == 0 else "")
        
        # Add value labels
        for i, (direct, indirect) in enumerate(valid_pairs):
            direct_sum = damage_totals[direct]
            indirect_sum = damage_totals[indirect]
            total = direct_sum + indirect_sum
            
            # Display direct values in middle of direct bar
            ax.text(i, direct_sum/2, f'{direct_sum:.0f}', ha='center', va='center', 
                   color='white', fontweight='bold')
            
            # Display indirect values in middle of indirect bar
            ax.text(i, direct_sum + indirect_sum/2, f'{indirect_sum:.0f}', ha='center', 
                   va='center', color='white', fontweight='bold')
            
            # Display total on top
            ax.text(i, total + 0.5, f'Total: {total:.0f}', ha='center', va='bottom')
        
        ax.set_xlabel('Damage Type', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Direct vs Indirect Damage Distribution', fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        
    plt.tight_layout()
    save_or_show_plot(fig, 'damage_distribution')
    return fig

def plot_gvalues(gvalue_results):
    """Create a horizontal bar plot showing G-values for chemical species.
    
    Args:
        gvalue_results: Dictionary with chemical species as keys, each containing
                       'value' and 'error' keys for the G-value and its uncertainty.
    """
    # Convert results to dataframe format
    data = {
        'Species': [],
        'GValue': [],
        'Error': []
    }
    
    for species, result in gvalue_results.items():
        if 'value' in result:
            data['Species'].append(species)
            data['GValue'].append(result['value'])
            data['Error'].append(result['error'] / 2)  # Convert from 2σ to 1σ for error bars
    
    # Convert to DataFrame and sort by GValue
    df = pd.DataFrame(data)
    df = df.sort_values('GValue')
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 8))  # Taller figure for better species label display
    
    # Use categorical colormap based on value
    colors = plt.cm.viridis(df['GValue'] / df['GValue'].max())
    
    bars = ax.barh(
        df['Species'],
        df['GValue'],
        xerr=df['Error'],
        align='center',
        ecolor='black',
        capsize=3,
        color=colors,
        alpha=0.7
    )
    
    # Add value annotations
    for i, (bar, value) in enumerate(zip(bars, df['GValue'])):
        ax.text(value + df['Error'][i] + 0.05, i, 
                f'{value:.3f}', 
                va='center', fontsize=9)
    
    ax.set_xlabel('G-Value (molecules / 100 eV)', fontsize=12)
    ax.set_title('Chemical Species Production (G-Values)', fontsize=14)
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add a colorbar legend
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                              norm=plt.Normalize(vmin=0, vmax=df['GValue'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.01)
    cbar.set_label('G-Value Magnitude', fontsize=10)
    
    plt.tight_layout()
    save_or_show_plot(fig, 'gvalues_distribution')
    return fig

def plot_chemical_species_violin(all_cell_results):
    """Create violin plot for chemical species G-values across cells.
    
    Args:
        all_cell_results: List of dictionaries, each containing results for one cell
                         with GValues dictionary containing species data.
    """
    # Prepare data
    species_data = defaultdict(list)
    for cell_results in all_cell_results:
        for species, data in cell_results['GValues'].items():
            species_data[species].append(data['value'])
    
    # Convert to DataFrame and create violin plot
    df = pd.DataFrame(species_data)
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Get a colormap for different species
    colors = cm.tab10(np.linspace(0, 1, len(df.columns)))
    
    violin_parts = ax.violinplot(
        [df[col].values for col in df.columns],
        showmeans=True, 
        showmedians=True,
        vert=True
    )
    
    # Customize violin plot
    ax.set_xticks(range(1, len(df.columns) + 1))
    ax.set_xticklabels(df.columns, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('G-Value (molecules/100eV)', fontsize=12)
    ax.set_title('Distribution of G-Values Across Cells', fontsize=14)
    
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    
    # Color the violins
    for i, (pc, color) in enumerate(zip(violin_parts['bodies'], colors)):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
    
    # Color the median and mean lines
    for partname, part in violin_parts.items():
        if partname != 'bodies':
            if partname == 'cmeans':
                part.set_edgecolor('red')
                part.set_linewidth(1.5)
            elif partname == 'cmedians':
                part.set_edgecolor('black')
                part.set_linewidth(1.5)
            
    # Add a legend
    legend_elements = [
        Line2D([0], [0], color='red', lw=1.5, label='Mean'),
        Line2D([0], [0], color='black', lw=1.5, label='Median')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Add annotations for mean values
    for i, col in enumerate(df.columns):
        mean_val = df[col].mean()
        ax.text(i+1, mean_val, f'{mean_val:.3f}', 
                ha='center', va='bottom', fontsize=8, 
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
    
    plt.tight_layout()
    save_or_show_plot(fig, 'chemical_species_violin')
    return fig

def plot_dna_damage_violin(all_cell_results):
    """Create violin plot for DNA damage across cells.
    
    Args:
        all_cell_results: List of dictionaries, each containing results for one cell
                         with DNADamage dictionary containing damage statistics.
    """
    # Define damage types and their attributes
    damage_types = {
        'DSB': {'color': '#1f77b4', 'label': 'Double Strand Breaks'},
        'SSB': {'color': '#ff7f0e', 'label': 'Single Strand Breaks'},
        'SB': {'color': '#2ca02c', 'label': 'Strand Breaks'},
        'BD': {'color': '#d62728', 'label': 'Base Damage'}
    }
    
    # Collect data
    damage_data = defaultdict(list)
    for cell_results in all_cell_results:
        for damage_type in damage_types.keys():
            damage_data[damage_type].append(cell_results['DNADamage'][damage_type])
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 7))
    df = pd.DataFrame(damage_data)
    
    # Create violin plots with enhanced styling
    violin_parts = ax.violinplot(
        [df[col].values for col in df.columns], 
        showmeans=True, 
        showmedians=True,
        vert=True
    )
    
    # Customize plot
    ax.set_xticks(range(1, len(df.columns) + 1))
    ax.set_xticklabels([damage_types[col]['label'] for col in df.columns], fontsize=11)
    ax.set_ylabel('Number of Events', fontsize=12)
    ax.set_title('Distribution of DNA Damage Events Across Cells', fontsize=14)
    
    # Color the violins
    for i, pc in enumerate(violin_parts['bodies']):
        damage_type = list(damage_types.keys())[i]
        pc.set_facecolor(damage_types[damage_type]['color'])
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
    
    # Style mean and median lines
    violin_parts['cmeans'].set_edgecolor('red')
    violin_parts['cmeans'].set_linewidth(1.5)
    violin_parts['cmedians'].set_edgecolor('black')
    violin_parts['cmedians'].set_linewidth(1.5)
    
    # Add a legend for mean and median
    legend_elements = [
        Line2D([0], [0], color='red', lw=1.5, label='Mean'),
        Line2D([0], [0], color='black', lw=1.5, label='Median')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Add mean value annotations
    for i, col in enumerate(df.columns):
        mean_val = df[col].mean()
        ax.text(i+1, mean_val, f'{mean_val:.1f}', 
                ha='center', va='bottom', fontsize=9,
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
    
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    plt.tight_layout()
    save_or_show_plot(fig, 'dna_damage_violin')
    return fig

def display_multicell_results(all_cell_results, multicell_stats, save_plots=False, output_path=None):
    """Display comprehensive results from multicell analysis including tables and plots.
    
    Args:
        all_cell_results: List of results from each cell
        multicell_stats: Aggregated statistics across cells
        save_plots: If True, saves plots to files instead of displaying them.
        output_path: Optional path where to save results. If None, uses current directory.
    """
    # Create DataFrame with results table
    data = []

    # Define all columns
    columns = [
        'Cell',
        'DoseToNucl_ph2 (Gy)', 
        'DoseToNucl_ph3 (Gy)',
        'Energy to Cell (MeV)',
        'NP electrons',
        'DSB', 'SSB', 'SB', 'BD'
    ]

    # Add G-Value columns for each species
    species_list = list(all_cell_results[0]['GValues'].keys())
    gvalue_columns = [f'G({species})' for species in species_list]
    columns.extend(gvalue_columns)

    # Collect data for each cell
    for i, cell_results in enumerate(all_cell_results):
        row = [
            f'{i+1}',
            cell_results['DoseToNucl_ph2']['value'],
            cell_results['DoseToNucl_ph3']['value'],
            cell_results['Ecell']['value'],
            cell_results['NP_el']['value'],
            cell_results['DNADamage']['DSB'],
            cell_results['DNADamage']['SSB'],
            cell_results['DNADamage']['SB'],
            cell_results['DNADamage']['BD']
        ]
        # Add G-Values
        for species in species_list:
            row.append(cell_results['GValues'][species]['value'])
        data.append(row)

    # Add mean values row
    mean_row = [
        'Mean',
        multicell_stats['DoseToNucl_ph2']['mean'],
        multicell_stats['DoseToNucl_ph3']['mean'],
        multicell_stats['Ecell']['mean'],
        multicell_stats['NP_el']['mean'],
        multicell_stats['DNADamage']['DSB']['mean'],
        multicell_stats['DNADamage']['SSB']['mean'],
        multicell_stats['DNADamage']['SB']['mean'],
        multicell_stats['DNADamage']['BD']['mean']
    ]
    # Add G-Value means
    for species in species_list:
        mean_row.append(multicell_stats['GValues'][species]['mean'])
    data.append(mean_row)

    # Add standard deviation row
    error_row = [
        'Uncertainty',
        multicell_stats['DoseToNucl_ph2']['error'],
        multicell_stats['DoseToNucl_ph3']['error'],
        multicell_stats['Ecell']['error'],
        multicell_stats['NP_el']['error'],
        multicell_stats['DNADamage']['DSB']['error'],
        multicell_stats['DNADamage']['SSB']['error'],
        multicell_stats['DNADamage']['SB']['error'],
        multicell_stats['DNADamage']['BD']['error']
    ]
    # Add G-Value standard deviations
    for species in species_list:
        error_row.append(multicell_stats['GValues'][species]['error'])
    data.append(error_row)

    # Create DataFrame and format display
    results_df = pd.DataFrame(data, columns=columns)

    # Format numbers with appropriate precision
    def format_value(x):
        if pd.isna(x):
            return ''
        if isinstance(x, (int, np.integer)) or (isinstance(x, float) and x.is_integer()):
            return f'{int(x)}'
        if isinstance(x, float):
            if abs(x) < 1e-4 or abs(x) > 1e4:
                return f'{x:.2e}'
            return f'{x:.4f}'
        return str(x)

    # Convert integer columns to int type and format all columns
    formatted_df = pd.DataFrame()
    for col in results_df.columns:
        if col in ['NP electrons', 'DSB', 'SSB', 'SB', 'BD']:
            formatted_df[col] = results_df[col].astype('float').round().astype('Int64')
        else:
            formatted_df[col] = results_df[col].apply(format_value)

    # Add separator lines for mean and uncertainty rows
    print("\nResults Table:")
    print("=" * 120)  # Adjust width as needed
    
    # Print the header
    print(formatted_df.columns.str.ljust(15).str.cat(sep=' '))
    print("-" * 120)  # Separator line
    
    # Print the data rows
    for idx, row in formatted_df.iterrows():
        if idx == len(formatted_df) - 2:  # Before mean row
            print("-" * 120)
        print(row.str.ljust(15).str.cat(sep=' '))
        if idx == len(formatted_df) - 1:  # After last row
            print("=" * 120)
    
    # Save table to CSV if in terminal mode
    if not is_jupyter():
        _, tables_dir = get_output_dirs(output_path)
        csv_filename = os.path.join(tables_dir, 'multicell_results.csv')
        formatted_df.to_csv(csv_filename, index=False)
        print(f"\nTable saved to {csv_filename}")

    # Plot distributions
    plot_chemical_species_violin(all_cell_results)
    plot_dna_damage_violin(all_cell_results)
    
    if not is_jupyter():
        print("\nPlots have been saved in the 'plots' directory:"
              "\n- chemical_species_violin.png"
              "\n- dna_damage_violin.png")

def create_enhancement_bar_plot(data_list, labels, errors_list, title, colors=None, scenario_labels=None):
    """Create a bar plot for a specific enhancement category with multiple scenarios.
    
    Args:
        data_list: List of lists, where each inner list contains enhancement ratio values for a scenario
        labels: List of labels for each bar
        errors_list: List of lists, where each inner list contains error values for a scenario
        title: Title for the plot
        colors: List of colors for each scenario (will use default colors if None)
        scenario_labels: List of labels for each scenario (will use "Scenario X" if None)
    """
    if not isinstance(data_list[0], list):
        # Handle the case of a single scenario (backward compatibility)
        data_list = [data_list]
        errors_list = [errors_list]
        if scenario_labels is None:
            scenario_labels = [""]
    
    n_scenarios = len(data_list)
    n_categories = len(labels)
    
    # Set up colors if not provided
    if colors is None:
        cmap = plt.cm.tab10
        colors = [cmap(i/10) for i in range(n_scenarios)]
    
    # Set up scenario labels if not provided
    if scenario_labels is None:
        scenario_labels = [f"Scenario {i+1}" for i in range(n_scenarios)]
    
    # Bar width based on number of scenarios
    bar_width = 0.7 / n_scenarios
    
    # Set up the figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Position adjustment for each scenario's bars
    positions = [np.arange(n_categories) - 0.35 + (i + 0.5) * bar_width for i in range(n_scenarios)]
    
    # Plot bars for each scenario
    bars_list = []
    for i in range(n_scenarios):
        data = data_list[i]
        errors = errors_list[i]
        pos = positions[i]
        
        # Skip scenarios with no data
        if len(data) == 0:
            continue
            
        bars = ax.bar(pos, data, bar_width, yerr=errors, capsize=3, 
                     color=colors[i], alpha=0.7, label=scenario_labels[i])
        bars_list.append(bars)
    
    # Add horizontal line at y=1
    ax.axhline(y=1, color='k', linestyle='--', alpha=0.3)
    
    # Customize plot
    ax.set_xticks(np.arange(n_categories))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Enhancement Ratio', fontsize=12)
    ax.set_title(f'{title} Enhancement', fontsize=14)
    
    # Add legend if multiple scenarios
    if n_scenarios > 1:
        ax.legend(fontsize=10)
    
    # Add value labels for each bar
    for i, bars in enumerate(bars_list):
        for j, (bar, v, err) in enumerate(zip(bars, data_list[i], errors_list[i])):
            y_pos = v + err + 0.05
            ax.text(bar.get_x() + bar.get_width()/2, y_pos, 
                   f'{v:.2f}', ha='center', va='bottom', fontsize=8,
                   rotation=45 if n_scenarios > 1 else 0)
    
    # Add grid
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    return fig

def extract_enhancement_data(enhancement_results, category):
    """Extract enhancement data for a specific category.
    
    Args:
        enhancement_results: The output from compute_enhancement_ratios
        category: One of 'dose_energy', 'gvalues', 'dna_damage', 'complexity',
                 'dna_damage_per_gy', or 'complexity_per_gy'
        
    Returns:
        Tuple of (data, labels, errors)
    """
    if category == 'dose_energy':
        # Extract dose and energy data
        data = []
        labels = []
        errors = []
        
        for key in ['DoseToNucl_ph2', 'DoseToNucl_ph3', 'Ecell']:
            if key in enhancement_results['simple_quantities']:
                result = enhancement_results['simple_quantities'][key]
                data.append(result['ratio'])
                errors.append(result['uncertainty'])
                
                # Use more descriptive labels
                display_name = {
                    'DoseToNucl_ph2': 'Dose to Nucleus\n(Phase 2)',
                    'DoseToNucl_ph3': 'Dose to Nucleus\n(Phase 3)',
                    'Ecell': 'Energy to Cell'
                }[key]
                labels.append(display_name)
        
    elif category == 'gvalues':
        # Extract G-Values data
        data = []
        labels = []
        errors = []
        
        for species in enhancement_results['GValues']:
            result = enhancement_results['GValues'][species]
            data.append(result['ratio'])
            errors.append(result['uncertainty'])
            labels.append(f'G({species})')
    
    elif category == 'dna_damage':
        # Extract DNA Damage data
        data = []
        labels = []
        errors = []
        
        for key in ['DSB', 'SSB', 'SB', 'BD']:
            if key in enhancement_results['DNADamage']:
                result = enhancement_results['DNADamage'][key]
                data.append(result['ratio'])
                errors.append(result['uncertainty'])
                
                # Use more descriptive labels
                display_name = {
                    'DSB': 'Double Strand\nBreaks',
                    'SSB': 'Single Strand\nBreaks',
                    'SB': 'Strand\nBreaks',
                    'BD': 'Base\nDamage'
                }[key]
                labels.append(display_name)
    
    elif category == 'complexity':
        # Extract Complexity data
        data = []
        labels = []
        errors = []
        complexity_pattern = re.compile(r'Complexity\d+')
        
        # Collect complexity data
        for damage_type in enhancement_results['DNADamage']:
            if complexity_pattern.match(damage_type):
                result = enhancement_results['DNADamage'][damage_type]
                if result['ratio'] is not None:
                    data.append(result['ratio'])
                    errors.append(result['uncertainty'])
                    complexity_number = re.search(r'\d+', damage_type).group()
                    labels.append(f'Complexity {complexity_number}')
        
        # Sort by complexity number
        if data:
            # Sort all lists by the complexity number
            sort_indices = sorted(range(len(labels)), 
                                key=lambda i: int(re.search(r'\d+', labels[i]).group()))
            
            data = [data[i] for i in sort_indices]
            labels = [labels[i] for i in sort_indices]
            errors = [errors[i] for i in sort_indices]
    
    elif category == 'dna_damage_per_gy':
        # Extract DNA Damage per Gy data
        data = []
        labels = []
        errors = []
        
        if 'DNADamage_per_Gy' in enhancement_results:
            for key in ['DSB', 'SSB', 'SB', 'BD']:
                if key in enhancement_results['DNADamage_per_Gy']:
                    result = enhancement_results['DNADamage_per_Gy'][key]
                    data.append(result['ratio'])
                    errors.append(result['uncertainty'])
                    
                    # Use more descriptive labels
                    display_name = {
                        'DSB': 'Double Strand\nBreaks/Gy',
                        'SSB': 'Single Strand\nBreaks/Gy',
                        'SB': 'Strand\nBreaks/Gy',
                        'BD': 'Base\nDamage/Gy'
                    }[key]
                    labels.append(display_name)
    
    elif category == 'complexity_per_gy':
        # Extract Complexity per Gy data
        data = []
        labels = []
        errors = []
        complexity_pattern = re.compile(r'Complexity\d+')
        
        if 'DNADamage_per_Gy' in enhancement_results:
            # Collect complexity data
            for damage_type in enhancement_results['DNADamage_per_Gy']:
                if complexity_pattern.match(damage_type):
                    result = enhancement_results['DNADamage_per_Gy'][damage_type]
                    if result['ratio'] is not None:
                        data.append(result['ratio'])
                        errors.append(result['uncertainty'])
                        complexity_number = re.search(r'\d+', damage_type).group()
                        labels.append(f'Complexity {complexity_number}/Gy')
            
            # Sort by complexity number
            if data:
                # Sort all lists by the complexity number
                sort_indices = sorted(range(len(labels)), 
                                    key=lambda i: int(re.search(r'\d+', labels[i]).group()))
                
                data = [data[i] for i in sort_indices]
                labels = [labels[i] for i in sort_indices]
                errors = [errors[i] for i in sort_indices]
    
    return data, labels, errors

def plot_multi_enhancement_categories(enhancement_results_list):
    """Create plots for all enhancement categories showing multiple scenarios.
    
    Args:
        enhancement_results_list: List of outputs from compute_enhancement_ratios
    """
    # Define categories and their colors
    categories = [
        ('dose_energy', 'Dose and Energy'),
        ('gvalues', 'G-Values'),
        ('dna_damage', 'DNA Damage'),
        ('complexity', 'Complexity'),
        ('dna_damage_per_gy', 'DNA Damage per Gy'),
        ('complexity_per_gy', 'Complexity per Gy')
    ]
    
    # Get scenario labels
    scenario_labels = [er.get('scenario_label', f'Scenario {i+1}') 
                      for i, er in enumerate(enhancement_results_list)]
    
    # Set colors using colormap
    cmap = plt.cm.tab10
    colors = [cmap(i/10) for i in range(len(enhancement_results_list))]
    
    # Process each category
    all_figures = []
    for category_key, title in categories:
        all_data = []
        all_errors = []
        all_labels = set()
        
        # Collect data from each scenario for the current category
        for er in enhancement_results_list:
            data, labels, errors = extract_enhancement_data(er, category_key)
            all_data.append(data)
            all_errors.append(errors)
            all_labels.update(labels)
        
        # If no data for this category, skip
        if all(len(data) == 0 for data in all_data):
            continue
            
        # Get common labels across all scenarios (to ensure consistent ordering)
        if category_key in ['complexity', 'complexity_per_gy']:
            # For complexity categories, sort by number
            common_labels = sorted(list(all_labels), 
                                 key=lambda x: int(re.search(r'\d+', x).group()))
        else:
            # For other categories, use the predefined order
            predefined_labels = {
                'dose_energy': ['Dose to Nucleus\n(Phase 2)', 'Dose to Nucleus\n(Phase 3)', 'Energy to Cell'],
                'dna_damage': ['Double Strand\nBreaks', 'Single Strand\nBreaks', 'Strand\nBreaks', 'Base\nDamage'],
                'dna_damage_per_gy': ['Double Strand\nBreaks/Gy', 'Single Strand\nBreaks/Gy', 
                                     'Strand\nBreaks/Gy', 'Base\nDamage/Gy'],
                'gvalues': sorted(list(all_labels))  # Sort alphabetically for G-values
            }
            common_labels = predefined_labels.get(category_key, sorted(list(all_labels)))
        
        # Reorder and align data for each scenario based on common labels
        aligned_data = []
        aligned_errors = []
        
        for i, (data, errors, labels_used) in enumerate(zip(all_data, all_errors, 
                                                 [extract_enhancement_data(er, category_key)[1] 
                                                  for er in enhancement_results_list])):
            scenario_data = []
            scenario_errors = []
            
            # For each common label, find corresponding data or use NaN
            for label in common_labels:
                if label in labels_used:
                    idx = labels_used.index(label)
                    scenario_data.append(data[idx])
                    scenario_errors.append(errors[idx])
                else:
                    scenario_data.append(float('nan'))
                    scenario_errors.append(0)
            
            aligned_data.append(scenario_data)
            aligned_errors.append(scenario_errors)
        
        # Create the multi-scenario bar plot for this category
        fig = create_enhancement_bar_plot(aligned_data, common_labels, aligned_errors, 
                                        title, colors=colors, scenario_labels=scenario_labels)
        all_figures.append(fig)
    
    return all_figures

def plot_all_enhancement_categories(enhancement_results):
    """Create plots for all enhancement categories for a single scenario.
    
    Args:
        enhancement_results: The output from compute_enhancement_ratios
    """
    # Define categories and their colors
    categories = [
        ('dose_energy', '#3498db', 'Dose and Energy'),            # blue
        ('gvalues', '#2ecc71', 'G-Values'),                       # green
        ('dna_damage', '#e74c3c', 'DNA Damage'),                  # red
        ('complexity', '#9b59b6', 'Complexity'),                  # purple
        ('dna_damage_per_gy', '#e67e22', 'DNA Damage per Gy'),    # orange
        ('complexity_per_gy', '#8e44ad', 'Complexity per Gy')     # dark purple
    ]
    
    # Get scenario label if available
    scenario = enhancement_results.get('scenario_label', '')
    
    # Create plots for each category
    all_figures = []
    for category, color, title in categories:
        data, labels, errors = extract_enhancement_data(enhancement_results, category)
        if data:  # Only create plot if we have data
            fig = create_enhancement_bar_plot(data, labels, errors, title, 
                                            colors=[color], 
                                            scenario_labels=[scenario] if scenario else None)
            all_figures.append(fig)
    
    return all_figures

def display_enhancement_table_grouped(enhancement_results):
    """Display enhancement ratios grouped in separate tables.
    
    Args:
        enhancement_results: The output from compute_enhancement_ratios
        
    Returns:
        Dictionary of DataFrames with different enhancement categories
    """
    # Define column names for all tables
    columns = ['Quantity', 'Enhancement Ratio', 'Uncertainty']
    
    # ---- Group 1: Dose and Energy ----
    data_dose_energy = []
    for key in ['DoseToNucl_ph2', 'DoseToNucl_ph3', 'Ecell']:
        if key in enhancement_results['simple_quantities']:
            display_name = {
                'DoseToNucl_ph2': 'Dose to Nucleus (Phase 2)',
                'DoseToNucl_ph3': 'Dose to Nucleus (Phase 3)',
                'Ecell': 'Energy to Cell'
            }[key]
            result = enhancement_results['simple_quantities'][key]
            data_dose_energy.append([
                display_name,
                result['ratio'],
                result['uncertainty']
            ])
    
    # ---- Group 2: G-Values ----
    data_gvalues = []
    for species in enhancement_results['GValues']:
        result = enhancement_results['GValues'][species]
        data_gvalues.append([
            f'G-Value ({species})',
            result['ratio'],
            result['uncertainty']
        ])
    
    # ---- Group 3: DNA Damage (main types) ----
    data_dna_damage = []
    dna_damage_types = {
        'DSB': 'Double Strand Breaks',
        'SSB': 'Single Strand Breaks',
        'SB': 'Strand Breaks',
        'BD': 'Base Damage'
    }
    
    for damage_type, display_name in dna_damage_types.items():
        if damage_type in enhancement_results['DNADamage']:
            result = enhancement_results['DNADamage'][damage_type]
            if result['ratio'] is not None:
                data_dna_damage.append([
                    display_name,
                    result['ratio'],
                    result['uncertainty']
                ])
    
    # ---- Group 4: Complexity ----
    data_complexity = []
    complexity_pattern = re.compile(r'Complexity\d+')
    
    for damage_type in enhancement_results['DNADamage']:
        if complexity_pattern.match(damage_type):
            result = enhancement_results['DNADamage'][damage_type]
            if result['ratio'] is not None:
                data_complexity.append([
                    f'{damage_type} Damage',
                    result['ratio'],
                    result['uncertainty']
                ])
    
    # ---- Group 5: DNA Damage per Gy ----
    data_dna_damage_per_gy = []
    if 'DNADamage_per_Gy' in enhancement_results:
        for damage_type, display_name in dna_damage_types.items():
            if damage_type in enhancement_results['DNADamage_per_Gy']:
                result = enhancement_results['DNADamage_per_Gy'][damage_type]
                if result['ratio'] is not None:
                    data_dna_damage_per_gy.append([
                        f'{display_name} per Gy',
                        result['ratio'],
                        result['uncertainty']
                    ])
    
    # ---- Group 6: Complexity per Gy ----
    data_complexity_per_gy = []
    if 'DNADamage_per_Gy' in enhancement_results:
        for damage_type in enhancement_results['DNADamage_per_Gy']:
            if complexity_pattern.match(damage_type):
                result = enhancement_results['DNADamage_per_Gy'][damage_type]
                if result['ratio'] is not None:
                    data_complexity_per_gy.append([
                        f'{damage_type} Damage per Gy',
                        result['ratio'],
                        result['uncertainty']
                    ])
    
    # Sort complexity by number
    data_complexity.sort(key=lambda x: int(re.search(r'\d+', x[0]).group()))
    if data_complexity_per_gy:
        data_complexity_per_gy.sort(key=lambda x: int(re.search(r'\d+', x[0]).group()))
    
    def format_and_style_df(data, title):
        if not data:  # Skip empty data
            return None
            
        df = pd.DataFrame(data, columns=columns)
        
        # Format the ratios and uncertainties
        def format_ratio(x):
            if isinstance(x, float):
                return f"{x:.3f}"
            return str(x)
        
        # Apply formatting to columns
        formatted_df = df.copy()
        formatted_df['Enhancement Ratio'] = formatted_df['Enhancement Ratio'].apply(format_ratio)
        formatted_df['Uncertainty'] = formatted_df['Uncertainty'].apply(format_ratio)
        
        # Save to CSV if in terminal mode
        if not is_jupyter():
            os.makedirs('tables', exist_ok=True)
            safe_title = title.lower().replace(" ", "_").replace("/", "_")
            csv_filename = f'tables/enhancement_{safe_title}.csv'
            df.to_csv(csv_filename, index=False)
            print(f"\nEnhancement table for {title} saved to {csv_filename}")
        
        return formatted_df
    
    # Create and format all DataFrames
    dfs = {
        "Dose and Energy Enhancement": format_and_style_df(data_dose_energy, "Dose and Energy Enhancement"),
        "G-Values Enhancement": format_and_style_df(data_gvalues, "G-Values Enhancement"),
        "DNA Damage Enhancement": format_and_style_df(data_dna_damage, "DNA Damage Enhancement"),
        "Complexity Enhancement": format_and_style_df(data_complexity, "Complexity Enhancement"),
        "DNA Damage per Gy Enhancement": format_and_style_df(data_dna_damage_per_gy, "DNA Damage per Gy Enhancement"),
        "Complexity per Gy Enhancement": format_and_style_df(data_complexity_per_gy, "Complexity per Gy Enhancement")
    }
    
    return {k: v for k, v in dfs.items() if v is not None}

def display_enhancement_table(enhancement_results):
    """Original function that displays all enhancement ratios in a single table.
    
    Args:
        enhancement_results: The output from compute_enhancement_ratios
        
    Returns:
        Tuple of (DataFrame, StyledDataFrame) with all enhancement ratios
    """
    data = []
    columns = ['Quantity', 'Enhancement Ratio', 'Uncertainty']
    
    # Physical quantities
    for key in ['DoseToNucl_ph2', 'DoseToNucl_ph3', 'Ecell']:
        if key in enhancement_results['simple_quantities']:
            display_name = {
                'DoseToNucl_ph2': 'Dose to Nucleus (Phase 2)',
                'DoseToNucl_ph3': 'Dose to Nucleus (Phase 3)',
                'Ecell': 'Energy to Cell'
            }[key]
            result = enhancement_results['simple_quantities'][key]
            data.append([
                display_name,
                result['ratio'],
                result['uncertainty']
            ])
    
    # DNA Damage
    damage_display = {
        'DSB': 'Double Strand Breaks',
        'SSB': 'Single Strand Breaks',
        'SB': 'Strand Breaks',
        'BD': 'Base Damage'
    }
    # Add complexity levels
    for i in range(2, 15):
        damage_display[f'Complexity{i}'] = f'Complexity {i}'
    
    for damage_type, display_name in damage_display.items():
        if damage_type in enhancement_results['DNADamage']:
            result = enhancement_results['DNADamage'][damage_type]
            if result['ratio'] is not None:
                data.append([
                    display_name,
                    result['ratio'],
                    result['uncertainty']
            ])
    
    # DNA Damage per Gy
    if 'DNADamage_per_Gy' in enhancement_results:
        for damage_type, display_name in damage_display.items():
            if damage_type in enhancement_results['DNADamage_per_Gy']:
                result = enhancement_results['DNADamage_per_Gy'][damage_type]
                if result['ratio'] is not None:
                    data.append([
                        f'{display_name} per Gy',
                        result['ratio'],
                        result['uncertainty']
                ])
    
    # G-Values
    for species in enhancement_results['GValues']:
        result = enhancement_results['GValues'][species]
        data.append([
            f'G-Value ({species})',
            result['ratio'],
            result['uncertainty']
        ])
    
    # Create DataFrame
    df = pd.DataFrame(data, columns=columns)
    
    # Format the ratios and uncertainties
    def format_ratio(x):
        if isinstance(x, float):
            return f"{x:.3f}"
        return str(x)
    
    styled_df = df.style.format({
        'Enhancement Ratio': format_ratio,
        'Uncertainty': format_ratio
    })
    
    # Add coloring based on enhancement (>1 is blue, <1 is red)
    def color_enhancement(val):
        try:
            val = float(val)
            if val > 1:
                return 'color: blue'
            elif val < 1:
                return 'color: red'
        except:
            pass
        return ''
    
    styled_df = styled_df.applymap(color_enhancement, subset=['Enhancement Ratio'])
    
    return df, styled_df
