import matplotlib.pyplot as plt
import pandas as pd
import os


def plot_timeofdeath_histogram(csvpath, nbin=6, causes=["MitoticCatastrophe", "Apoptosis", "Necrosis"]):
    """
    Reads the cell death CSV and plots a histogram of time of death for each cause.
    Args:
        csvpath (str): Path to the CSV file.
        nbin (float): Bin width in hours.
        causes (list): List of death causes to plot.
    """
    # Read CSV
    df = pd.read_csv(csvpath)
    # Filter out rows with N/A or missing time of death
    df = df[df['Time of Death (h)'].apply(lambda x: str(x).replace('.', '', 1).isdigit())]
    df['Time of Death (h)'] = df['Time of Death (h)'].astype(float)
    # Prepare data for each cause
    data = []
    for cause in causes:
        times = df[df['Type of Death'] == cause]['Time of Death (h)'].values
        data.append(times)
    # Determine bin edges
    all_times = df['Time of Death (h)'].values
    if len(all_times) == 0:
        print("No valid time of death data to plot.")
        return
    min_time = 0
    max_time = max(all_times)
    bins = list(range(int(min_time), int(max_time) + int(nbin), int(nbin)))
    # Plot
    plt.figure(figsize=(10, 6))
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]  # Professional color palette
    nseries = len(data)
    # Use step histtype for clarity, with thicker lines
    for i, (series, cause) in enumerate(zip(data, causes)):
        plt.hist(series, bins=bins, histtype='step', linewidth=2.2, label=cause, color=colors[i % len(colors)])
        plt.hist(series, bins=bins, alpha=0.25, label=None, color=colors[i % len(colors)], stacked=True)
    plt.xlabel('Time of death [h]', fontsize=16)
    plt.ylabel('Number of cell deaths', fontsize=16)
    # plt.title('Distribution of Cell Death Times by Cause', fontsize=16, weight='bold', pad=15)
    plt.legend(fontsize=14, frameon=True, edgecolor='black')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout(pad=2)
    # Ensure the plot directory exists in the same directory as the CSV
    csv_dir = os.path.dirname(csvpath)
    os.makedirs(csv_dir, exist_ok=True)
    plot_path = os.path.join(csv_dir, 'time_of_death_histogram.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_DSBremaining(dsb_csv_path, std_type='std'):
    """
    Plot mean and std of DSB remaining over time from a CSV file.
    Args:
        dsb_csv_path (str): Path to the CSV file with columns 'Time (h)', 'Mean Fraction DSB Remaining',
            'Std Fraction DSB Remaining', and optionally 'Std_mean Fraction DSB Remaining'.
        std_type (str): 'std' for cell-to-cell std, 'std_mean' for std of the mean across repeats.
    """
    df = pd.read_csv(dsb_csv_path)
    time = df['Time (h)']
    mean = df['Mean Fraction DSB Remaining']
    if std_type == 'std_mean' and 'Std_mean Fraction DSB Remaining' in df.columns:
        std = df['Std_mean Fraction DSB Remaining']
        label = '±1 Std. Dev. (mean across repeats)'
        color = 'orange'
    else:
        std = df['Std Fraction DSB Remaining']
        label = '±1 Std. Dev. (cell-to-cell)'
        color = 'skyblue'
    plt.figure(figsize=(8,5))
    plt.plot(time, mean, color='navy', label='Average DSB Remaining')
    plt.fill_between(time, mean-std, mean+std, color=color, alpha=0.4, label=label)
    plt.xlabel('Time [h]', fontsize=16)
    plt.ylabel('Remaining DSB per cell', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
   # plt.title('DSB Remaining Over Time', fontsize=16, weight='bold')
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    # Ensure the plot directory exists in the same directory as the CSV
    csv_dir = os.path.dirname(dsb_csv_path)
    os.makedirs(csv_dir, exist_ok=True)
    plot_path = os.path.join(csv_dir, 'RemainingDSB.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_everything_in_dir(results_dir):
    """
    Plot all relevant CSV files in a directory.
    DSBremaining searches for file containing 'DSBremaining'
    timeofdeath_histogram searches for 'timeofdeath' in the filename.
    Args:
        results_dir (str): Directory containing the CSV files.
    """
    for filename in os.listdir(results_dir):
        if filename.endswith('.csv'):
            filepath = os.path.join(results_dir, filename)
            if 'DSBremaining' in filename:
                print(f"Plotting DSB remaining from {filename}")
                plot_DSBremaining(filepath, std_type='std')
            elif 'timeofdeath' in filename:
                print(f"Plotting time of death histogram from {filename}")
                plot_timeofdeath_histogram(filepath, nbin=6, causes=["MitoticCatastrophe", "Apoptosis", "Necrosis"])

results_dir = '/home/radiofisica/hector/mytopassimulations/TOPAS_CellsNPs/work/CellColony-med0-cell0/chronorepair'
plot_everything_in_dir(results_dir)
