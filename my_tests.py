import dnadamage_phsp_manager 
import sddparser
import pprint, os
import matplotlib.pyplot as plt


def test_dnadamagephsp_read():
    """
    Test the read_dnadamage_phase_space function.
    """
    base_path = "/home/radiofisica/hector/mytopassimulations/tests/run1-med1-cell1/DNADamage"

    # Read the phase space data.
    df = dnadamage_phsp_manager.read_dnadamage_phase_space(base_path)
    
    # Print the first few rows of the dataframe.
    print(df.head())
    
    # Calculate total damage statistics.
    dnadamage_phsp_manager.total_damage(df)



def test_dnadamagephsp_merge():
    """
    Test the merge_dnadamage_files function.
    """
    # Example usage:
    filebases = [
        "/home/radiofisica/hector/mytopassimulations/tests/run1-med1-cell1/DNADamage",
        "/home/radiofisica/hector/mytopassimulations/tests/run2-med1-cell1/DNADamage",
        # ... add more filebases as needed ...
    ]
    
    merged_data = dnadamage_phsp_manager.merge_dnadamage_files(filebases)
    
    # Print the first few rows of the merged data.
    print(merged_data.head())


def test_parseSDDFile():
    """
    Test the parseSDDFile function.
    """
    fileName = "/home/radiofisica/hector/mytopassimulations/tests/run1-med1-cell1/DNADamage_sdd.txt"
    verbose = False

    header, events = sddparser.parseSDDFileFlat(fileName, verbose)
    
    pp = pprint.PrettyPrinter(indent=2)
    print("=== Header ===")
    pp.pprint(header)
    print('\n')
    print("\n=== First 3 Events ===")
    for idx, event in enumerate(events[:3]):
        print(f"\nEvent #{idx + 1}:")
        pp.pprint(event)
        print('\n')

def test_multirun():
    from analize_cell_sim_results import multirun_processing
    from display_cell_sim_results import display_results, plot_damage_distribution, plot_gvalues

    # Set parameters for multirun processing
    nruns = 100
    filebase = '../TOPAS_CellsNPs/work/CellColony-med0-cell0/cell17'

    # Process all runs and get results
    Cell_results = multirun_processing(nruns, filebase)
    print("Cell results:")
    pprint.pprint(Cell_results)
    # Display results using imported functions
    display_results(Cell_results)

    # Plot damage distribution if DNA damage data is available
    if 'DNADamage' in Cell_results:
        plot_damage_distribution(Cell_results['DNADamage'])

    # Plot G-values if chemical species data is available    
   # if any('value' in data for data in Cell_results['GValues'].values()):
   #     plot_gvalues(Cell_results['GValues'])

     # Plot number of molecules if chemical species data is available    
   # if any('value' in data for data in Cell_results['NumberOfMolecules'].values()):
   #     plot_gvalues(Cell_results['NumberOfMolecules'])
    


def test_multicell_analysis():
    # Process conditions without nanoparticles
    from analize_cell_sim_results import read_multicell_json, multicell_processing, process_multicell_results, compute_enhancement_ratios
    from display_cell_sim_results import (
        display_results, 
        plot_damage_distribution,
        plot_gvalues, 
        display_multicell_results,
        plot_all_enhancement_categories,
        plot_multi_enhancement_categories,
        display_enhancement_table_grouped
    )

    # Set parameters for multicell processing
    n_cells = 40 # Number of cells to process
    n_runs = 100  # Number of runs per cell
    base_dir = '../TOPAS_CellsNPs/work/CellColony-med0-cell0'  # Base directory containing cell directories

    # Process all cells and their runs
   # all_cell_results_med_cell = multicell_processing(n_cells, n_runs, base_dir, save_json=True)

    # Read from the JSON file if it exists
    json_path = os.path.join(base_dir, 'multicell_results.json')
    all_cell_results_med_cell = read_multicell_json(json_path)


    # Compute statistics across cells
    multicell_stats_med_cell = process_multicell_results(all_cell_results_med_cell)

    # Display results for condition without nanoparticles
    display_multicell_results(all_cell_results_med_cell, multicell_stats_med_cell, output_path=base_dir)


def test_enhancement_ratios():
    from analize_cell_sim_results import multicell_processing, process_multicell_results, compute_enhancement_ratios, read_multicell_json
    from display_cell_sim_results import plot_all_enhancement_categories, plot_multi_enhancement_categories, display_enhancement_table_grouped, plot_multicell_categories

    list_multicell_stats = []

    n_cells = 40  # Number of cells to process
    n_runs = 100  # Number of runs per cell
    # Process conditions without nanoparticles  
    base_dir = '../TOPAS_CellsNPs/work/CellColony-med0-cell0'  # Base directory containing cell directories
    #all_cell_results_med0_cell0 = multicell_processing(n_cells, n_runs, base_dir)
    json_path = os.path.join(base_dir, 'multicell_results.json')
    all_cell_results_med0_cell0 = read_multicell_json(json_path)
    multicell_stats_med0_cell0 = process_multicell_results(all_cell_results_med0_cell0)
    list_multicell_stats.append(multicell_stats_med0_cell0)

    # Process conditions with nanoparticles 1mg/ml 
    base_dir = '../TOPAS_CellsNPs/work/CellColony-med1-cell1'  # Base directory containing cell directories
    #all_cell_results_med1_cell1 = multicell_processing(n_cells, n_runs, base_dir)
    json_path = os.path.join(base_dir, 'multicell_results.json')
    all_cell_results_med1_cell1 = read_multicell_json(json_path)
    multicell_stats_med1_cell1 = process_multicell_results(all_cell_results_med1_cell1)
    list_multicell_stats.append(multicell_stats_med1_cell1)

    # Process conditions with nanoparticles 5mg/ml 
    base_dir = '../TOPAS_CellsNPs/work/CellColony-med5-cell5'  # Base directory containing cell directories
    #all_cell_results_med1_cell1 = multicell_processing(n_cells, n_runs, base_dir)
    json_path = os.path.join(base_dir, 'multicell_results.json')
    all_cell_results_med5_cell5 = read_multicell_json(json_path)
    multicell_stats_med5_cell5 = process_multicell_results(all_cell_results_med5_cell5)
    list_multicell_stats.append(multicell_stats_med5_cell5)


    # Store the enhancement results in a list for comparison
    all_enhancement_results = []

    # Compute enhancement ratios for 1mg/ml NPs vs Control
    enhancement_1mg = compute_enhancement_ratios(
        multicell_stats_med1_cell1, 
        multicell_stats_med0_cell0,
        scenario_label="1mg/ml NPs"
    )
    all_enhancement_results.append(enhancement_1mg)

    # Compute enhancement ratios for 5mg/ml NPs vs Control
    enhancement_5mg = compute_enhancement_ratios(
        multicell_stats_med5_cell5, 
        multicell_stats_med0_cell0,
        scenario_label="5mg/ml NPs"
    )
    all_enhancement_results.append(enhancement_5mg)

    # Plot all categories with both scenarios in the same plots
    print("\nMulti-Scenario Enhancement Comparison")
    print("------------------------------------")
    for scenario in all_enhancement_results:
        lable = scenario.get('scenario_label', 'Unknown Scenario')
        print(f"Scenario: {lable}")
        #pprint.pprint(scenario)
        print("-" * (len("Enhancement Ratios: ") + len(scenario['scenario_label'])))
        dfs_dic = display_enhancement_table_grouped(scenario)
        for k,v in dfs_dic.items():
            print(f'\n{k}:')
            print(v)
    
    figures = plot_multi_enhancement_categories(all_enhancement_results)
    figures.append(plot_multicell_categories(list_multicell_stats,['0 mg/ml', '1 mg/ml', '5 mg/ml'] ))

    # Visualize enhancement ratios for 1mg/ml NP concentration vs control
    #scenario_label = enhancement_results.get('scenario_label')
    #print(f"\nVisualization of Enhancement Ratios: {scenario_label}")
    #print("-" * (len("Visualization of Enhancement Ratios: ") + len(scenario_label)))
    # Can use either the old function or the new one with a single scenario
    # Using new function that supports multiple scenarios
    #figures = plot_multi_enhancement_categories([enhancement_results])
    
    # Display all returned figures
    for fig in figures:
        # plt.figure(fig)
        plt.show()

    

# To run the test:
#test_parseSDDFile()
# test_multirun()
test_multicell_analysis()
# test_enhancement_ratios()