import os
from files_and_directory_manager import get_outputfile_paths
import numpy as np


def merge_combine_csv(output_file_paths, output_path):

    os.makedirs(output_path, exist_ok=True)

    # Merge EnergyDeposit or DoseToMedium files from different runs
    combined_sum = 0
    combined_count_in_bin = 0
    combined_second_moment = 0
    combined_histories_with_scorer_active = 0
    cont = 0
    var = 0.0
    lines_list = []
    for file_path in output_file_paths:
        # subfolder_name = f'run{run_number}'

        # Read EnergyDepositToNucleus.csv file
        # file_path = Path(f'{folder_path}/{subfolder_name}/{filename}').absolute()
        path = os.path.dirname(file_path)
        run_number = path.split('run')[-1]
        filename = os.path.basename(file_path)

        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Extract the data from the last line
        if len(lines) > 0:
            data_line = lines[-1]
            lines_list.append(f'{run_number} {data_line}')
            data_values = [float(value) for value in data_line.split(', ')]

            sum_dose, mean_dose, count_in_bin, second_moment, variance, std_dev, histories_with_scorer_active = data_values

            # Combine the output values
            combined_sum += sum_dose
            combined_count_in_bin += count_in_bin
            combined_second_moment += second_moment
            combined_histories_with_scorer_active += histories_with_scorer_active
            var = var + sum_dose * sum_dose
            cont = cont + 1
    # Calculate the combined mean, variance, and standard deviation
    combined_mean = combined_sum / combined_histories_with_scorer_active
    combined_variance = combined_second_moment / combined_histories_with_scorer_active - combined_mean ** 2
    combined_std_dev = np.sqrt(combined_variance)
    var = np.sqrt(var / cont - combined_sum * combined_sum / (cont * cont))

    # Write the combined results to a new output file
    with open(os.path.join(output_path, f'combined_{filename}'), "w") as f:
        for line in lines[:-1]:
            f.write(line)

        f.write(f"{combined_sum}, {combined_mean}, {combined_count_in_bin}, {combined_second_moment}, "
                f"{combined_variance}, {combined_std_dev}, {combined_histories_with_scorer_active}\n")

    # Write a results file with all the results from each job
    lines_list.sort(key=lambda x: int(x.split()[0]))
    with open(os.path.join(output_path, f'AllJobs_{filename}'), "w") as f:
        for line in lines[:-1]:
            f.write(line)
        for line in lines_list:
            f.write(line)

    unc_2sigma = combined_std_dev / np.sqrt(combined_histories_with_scorer_active)
    print('')
    print(f'File {filename:}')
    print(f'Number of results merged: {cont} out of {len(output_file_paths)}')
    print(f'Dose (Gy/hist): {combined_mean} +/- {unc_2sigma}')
    print(f'Sum Dose = {combined_sum / cont} +/- {var}')
    print('')

