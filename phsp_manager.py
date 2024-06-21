import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from mpl_toolkits.mplot3d import Axes3D


def generate_weighted_energy_histogram(folder, phsp_file, bins, min_energy, max_energy,
                                       particle="gammas", show=True, save_path=None, output_txt=None):
    # Initialize histogram counts and uncertainties
    hist_counts = np.zeros(bins)
    hist_uncertainties = np.zeros(bins)

    # Calculate bin width
    bin_width = (max_energy - min_energy) / bins

    # Particle type
    particle_type = 22  # by default photons
    if particle == "electrons":
        particle_type = 11

    # Update histogram counts using weights
    with open(os.path.join(folder, phsp_file + '.phsp'), 'r') as f:
        for line in f:
            cols = line.split()
            if int(cols[7]) == particle_type:
                energy = float(cols[5])  # Energy is at index 5
                weight = float(cols[6])  # Weight is at index 6
                bin_index = int((energy - min_energy) / bin_width)
                if bin_index == bins:  # Include the max energy value in the last bin
                    bin_index -= 1
                hist_counts[bin_index] += weight

    # Calculate uncertainties
    hist_uncertainties = np.sqrt(hist_counts)

    # Plot the histogram with weights
    bin_edges = np.linspace(min_energy, max_energy, bins + 1)
    plt.hist(bin_edges[:-1], bin_edges, weights=hist_counts)
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Weighted Counts')
    plt.title('Weighted Energy Histogram')

    if save_path:
        plt.savefig(os.path.join(folder, f'EnergyHistogram_{phsp_file}_{particle}.png'))

    if show:
        plt.show()

    plt.close()

    if output_txt:
        bin_centers = np.linspace(min_energy + bin_width / 2, max_energy - bin_width / 2, bins)
        with open(os.path.join(folder, f'EnergyHistogram_{phsp_file}_{particle}.txt'), 'w') as f:
            f.write('# Energy[MeV] Weighted_Counts Uncertainty\n')
            for i in range(bins):
                f.write(f'{bin_centers[i]} {hist_counts[i]} {hist_uncertainties[i]}\n')


def read_phsp_to_dataframe(phsp_file):
    # Example usage
    # phsp_file = 'example.phsp'
    # df = read_phsp_to_dataframe(phsp_file)
    # print(df)

    # Define the column names
    column_names = ['pos_x', 'pos_y', 'pos_z', 'dircos_x', 'dircos_y', 'energy', 'weight', 'particle_type',
                    'flag_neg_dircos', 'flag_first']

    # Read the whole phsp file and create a DataFrame
    df = pd.read_csv(phsp_file + '.phsp', sep='\s+', header=None, names=column_names)

    return df


def merge_phsp_files(phsp_files, output_path):
    if len(phsp_files) > 0:
        # Remove the output file if it exists
        filename = os.path.basename(phsp_files[0])
        output_phsp_file = os.path.join(output_path, filename)
        if os.path.exists(output_phsp_file):
            os.remove(output_phsp_file)
        with open(output_phsp_file, 'a') as outfile:
            for phsp_file in phsp_files:
                # Read and process the .phsp and .header files as needed
                # For example, read the .phsp file line by line, and extract information from the .header file
                # pass
                # Append the result_phsp.phsp file to the combined output file
                with open(phsp_file, 'r') as infile:
                    outfile.writelines(infile.readlines())


def merge_header_files(header_files, output_path):
    def read_val(str_line):
        val = str_line.split(":")[-1].strip().split(" ", 1)[0]
        val = float(val)

        return val

    kw_orig_hist = "Number of Original Histories:"
    kw_hist_phsp = "Number of Original Histories that Reached Phase Space"
    kw_part = "Number of Scored Particles"
    kw_e = "Number of e-"
    kw_gammas = "Number of gamma"
    kw_min_E_e = "Minimum Kinetic Energy of e-"
    kw_min_E_g = "Minimum Kinetic Energy of gamma"
    kw_max_E_e = "Maximum Kinetic Energy of e-"
    kw_max_E_g = "Maximum Kinetic Energy of gamma"

    num_orig_hist = 0
    num_hist_phsp = 0
    num_part = 0
    num_e = 0
    num_gammas = 0
    min_E_e = 1.0e30
    min_E_g = 1.0e30
    max_E_e = -1.0
    max_E_g = -1.0

    # Remove the output file if it exists
    filename = os.path.basename(header_files[0])
    output_file = os.path.join(output_path, filename)
    if os.path.exists(output_file):
        os.remove(output_file)

    with open(header_files[0], 'r') as file:
        combined_header = file.readlines()

    for header_file in header_files:
        with open(header_file, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if line.startswith(kw_orig_hist):
                    num_orig_hist += int(read_val(line))
                elif line.startswith(kw_hist_phsp):
                    num_hist_phsp += int(read_val(line))
                elif line.startswith(kw_part):
                    num_part += int(read_val(line))
                elif line.startswith(kw_e):
                    num_e += int(read_val(line))
                elif line.startswith(kw_gammas):
                    num_gammas += int(read_val(line))
                elif line.startswith(kw_min_E_e):
                    if read_val(line) < min_E_e:
                        min_E_e = read_val(line)
                elif line.startswith(kw_min_E_g):
                    if read_val(line) < min_E_g:
                        min_E_g = read_val(line)
                elif line.startswith(kw_max_E_e):
                    if read_val(line) > max_E_e:
                        max_E_e = read_val(line)
                elif line.startswith(kw_max_E_g):
                    if read_val(line) > max_E_g:
                        max_E_g = read_val(line)
    new_lines = []
    for line in combined_header:
        if line.startswith(kw_orig_hist):
            nline = f'{kw_orig_hist} {num_orig_hist}\n'
        elif line.startswith(kw_hist_phsp):
            nline = f'{kw_hist_phsp}: {num_hist_phsp}\n'
        elif line.startswith(kw_part):
            nline = f'{kw_part}: {num_part}\n'
        elif line.startswith(kw_e):
            nline = f'{kw_e}: {num_e}\n'
        elif line.startswith(kw_gammas):
            nline = f'{kw_gammas}: {num_gammas}\n'
        elif line.startswith(kw_min_E_e):
            nline = f'{kw_min_E_e}: {min_E_e} MeV\n'
        elif line.startswith(kw_min_E_g):
            nline = f'{kw_min_E_g}: {min_E_g} MeV\n'
        elif line.startswith(kw_max_E_e):
            nline = f'{kw_max_E_e}: {max_E_e} MeV\n'
        elif line.startswith(kw_max_E_g):
            nline = f'{kw_max_E_g}: {max_E_g} MeV\n'
        else:
            nline = line

        new_lines.append(nline)

    with open(output_file, 'w') as output_file:
        for line in new_lines:
            output_file.write(line)


def split_phsp_file(folder, input_file, n):
    # creo que está perdiendo algunas particulas, seguramente por la condicion de que cada chunk empiece con historia original
    if not os.path.exists(folder):
        os.makedirs(output_folder)

    with open(os.path.join(folder, f'{input_file}.phsp')) as infile:
        lines = infile.readlines()

    num_lines = len(lines)
    part_size = num_lines // n
    end_index = 0
    for i in range(n):
        start_index = end_index
        end_index = (i + 1) * part_size if i < n - 1 else num_lines

        # Ensure that the first particle of each part has a value of 1 for the last column
        while end_index < num_lines - 1 and int(lines[end_index + 1].split()[-1]) != 1:
            end_index += 1

        output_file = os.path.join(folder, f'work/run{i + 1}', f'{input_file}.phsp')
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        with open(output_file, 'w') as outfile:
            outfile.writelines(lines[start_index:end_index])

        generate_header(os.path.join(folder, f'work/run{i + 1}'), input_file)


def generate_header(folder, phsp_file):
    # Initialize variables
    total_histories = 0
    total_reached = 0
    total_scored = 0

    min_energy_electron = float('inf')
    min_energy_gamma = float('inf')
    max_energy_electron = 0
    max_energy_gamma = 0

    electron_count = 0
    gamma_count = 0

    # Particle type in PDG Format
    electron_pdg = 11
    gamma_pdg = 22

    # Read phsp file line by line
    with open(os.path.join(folder, phsp_file + '.phsp'), 'r') as f:
        for line in f:
            cols = line.split()
            x, y, z, dircos_x, dircos_y, energy, weight, particle_type, flag_neg_dircos, flag_first = cols
            energy = float(energy)
            particle_type = int(particle_type)
            flag_first = int(flag_first)

            total_scored += 1

            if flag_first:
                total_reached += 1

            if particle_type == electron_pdg:
                electron_count += 1
                min_energy_electron = min(min_energy_electron, energy)
                max_energy_electron = max(max_energy_electron, energy)
            elif particle_type == gamma_pdg:
                gamma_count += 1
                min_energy_gamma = min(min_energy_gamma, energy)
                max_energy_gamma = max(max_energy_gamma, energy)

    # Write the header file
    with open(os.path.join(folder, phsp_file + '.header'), 'w') as f:
        f.write("TOPAS ASCII Phase Space\n\n")
        f.write(f"Number of Original Histories: {total_histories}\n")
        f.write(f"Number of Original Histories that Reached Phase Space: {total_reached}\n")
        f.write(f"Number of Scored Particles: {total_scored}\n\n")
        f.write("Columns of data are as follows:\n")
        f.write(" 1: Position X [cm]\n 2: Position Y [cm]\n 3: Position Z [cm]\n 4: Direction Cosine X\n")
        f.write(" 5: Direction Cosine Y\n 6: Energy [MeV]\n 7: Weight\n 8: Particle Type (in PDG Format)\n")
        f.write(" 9: Flag to tell if Third Direction Cosine is Negative (1 means true)\n")
        f.write("10: Flag to tell if this is the First Scored Particle from this History (1 means true)\n\n")
        f.write(f"Number of e-: {electron_count}\n")
        f.write(f"Number of gamma: {gamma_count}\n\n")
        f.write(f"Minimum Kinetic Energy of e-: {min_energy_electron} MeV\n")
        f.write(f"Minimum Kinetic Energy of gamma: {min_energy_gamma} MeV\n\n")
        f.write(f"Maximum Kinetic Energy of e-: {max_energy_electron} MeV\n")
        f.write(f"Maximum Kinetic Energy of gamma: {max_energy_gamma} MeV\n")


def filter_phsp_by_particle(folder, phsp_file, particle):
    # particle type
    particle_type = 22  # by default fotons
    if particle == "electrons":
        particle_type = 11
    # Update histogram counts using weights
    with open(os.path.join(folder, phsp_file + '.phsp'), 'r') as f_source:
        with open(os.path.join(folder, phsp_file + '_' + particle + '.phsp'), 'w') as f_target:
            for line in f_source:
                cols = line.split()
                if int(cols[7]) == particle_type:
                    f_target.writelines(line)
    generate_header(folder, phsp_file + '_' + particle)


def plot_phsp_particle_positions(data, filename=None):
    # Extract position columns and energy
    positions = data[['pos_x', 'pos_y', 'pos_z']]
    energy = data['energy']

    # Create 3D scatter plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Normalize the energy values to the range [0, 1]
    energy_norm = (energy - energy.min()) / (energy.max() - energy.min())

    # Create a scatter plot with colormap based on energy
    scale_factor = 0.001  # for distance units conversion
    positions = scale_factor*positions
    scatter = ax.scatter(positions['pos_x'], positions['pos_y'],
                         positions['pos_z'], c=1000*energy, alpha=0.6, edgecolors='w', s=50, cmap='viridis')

    # Add a colorbar to show the energy values
    cbar = fig.colorbar(scatter, ax=ax, label='Energy [keV]')

    ax.set_xlabel('X [um]')
    ax.set_ylabel('Y [um]')
    ax.set_zlabel('Z [um]')
    ax.set_title('Particle Positions on Spherical Surface')

    # Save the figure to a file if the filename is provided
    if filename:
        plt.savefig(filename, dpi=300)
