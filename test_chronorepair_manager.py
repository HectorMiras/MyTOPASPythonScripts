# Import the module
from chronorepair_manager import setup_repair_simulation, run_repair_simulation, display_repair_results, plot_dsb_repair_kinetics

# Set up simulation parameters
sim_params = setup_repair_simulation(
    exposure_time=23,              # Hours
    simulation_time=48,            # Hours (usually matches exposure_time)
    time_steps=48,                 # Number of time steps for simulation
    nucleus_max_radius=4.65,       # Microns
    diffusion_model='free',        # Model for DNA fragment diffusion
    dsb_model='standard',          # Model for DSB repair
    ssb_model='standard',          # Model for SSB repair
    bd_model='standard',           # Model for base damage repair
    dose_rate_function='exponential',  # Function to model dose rate over time
    initial_dose_rate=0.13803/3600,    # Gy/hour (here converted to Gy/second)
    half_life=(59.39*24)*3600      # Half-life in seconds (here 59.39 days)
)

# Run the repair simulation
damage_path = '/home/hector/mytopassimulations/MGHsimulations/TOPAS_CellsNPs/work/only_results_CellColony-med5-cell5/' # Path to directory with cell/run structure
repair_results = run_repair_simulation(
    sim_params=sim_params,
    damage_path=damage_path,
    n_cells=10
)

# Display the results
display_repair_results(repair_results)

# Create visualization of repair kinetics
plot_dsb_repair_kinetics(repair_results, title='DSB Repair Kinetics for I-125 Exposure')