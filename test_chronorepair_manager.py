# Import the module
from chronorepair_manager import setup_repair_simulation, run_repair_simulation, display_repair_results, plot_dsb_repair_kinetics

# Set up simulation parameters
sim_params = setup_repair_simulation(
    exposure_time=23,
    nucleus_max_radius=4.65,
    initial_dose_rate=0.13803/3600,
    half_life=(59.39*24)*3600
)

# Run the repair simulation
damage_path = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/only_results_CellColony-med1-cell1/' # Path to directory with cell/run structure
repair_results = run_repair_simulation(
    sim_params=sim_params,
    damage_path=damage_path,
    n_cells=10
)

# Display the results
display_repair_results(repair_results)

# Create visualization of repair kinetics
plot_dsb_repair_kinetics(repair_results, title='DSB Repair Kinetics for I-125 Exposure')