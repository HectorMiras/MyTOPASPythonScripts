import os
from files_and_directory_manager import get_outputfile_paths
import numpy as np
import pandas as pd

# Function to merge the dose to medium csv files from parallel runs of the TOPAS_CellsNP simulations
def merge_CellsNP_csv(output_file_paths, output_path, append=False):

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
            if append:
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
    if append:
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


# Function to merge topas csv dose to medium outputs from parallel runs
def get_header_info(file_path):
    columns = []
    header_lines = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("#"):
            header_lines.append(line)
            if line.startswith("# DoseToMedium"):
                # Split the line by ":" and take the second element, which is the one with the column names
                column_part = line.split(":")[1]

                # Now split by spaces and use a list comprehension to strip each column name
                columns = [col_name.strip() for col_name in column_part.split() if col_name]

    # Get the numeric column elements from the last line
    data_cols = [data.strip() for data in lines[-1].split() if data]
    dif_len = len(data_cols) - len(columns)
    # Add coordinate index column names to the begining of the column names list
    if dif_len > 0:
        for ind in reversed(range(dif_len)):
            columns = [f"ind{ind}", *columns]

    return header_lines, columns


def merge_csv(output_file_paths, output_path):

    os.makedirs(output_path, exist_ok=True)

    # Merge EnergyDeposit or DoseToMedium files from different runs
    combined_sum = 0
    combined_count_in_bin = 0
    combined_second_moment = 0
    combined_histories_with_scorer_active = 0
    cont = 0
    var = 0.0
    lines_list = []

    header_lines, columns = get_header_info(output_file_paths[0])

    # Initialize an empty DataFrame to hold all the merged data
    merged_data = pd.DataFrame()
    cont = 0

    for file_path in output_file_paths:
        # subfolder_name = f'run{run_number}'

        # Read EnergyDepositToNucleus.csv file
        # file_path = Path(f'{folder_path}/{subfolder_name}/{filename}').absolute()
        path = os.path.dirname(file_path)
        run_number = path.split('run')[-1]
        filename = os.path.basename(file_path)

        with open(file_path, 'r') as f:
            lines = f.readlines()

        current_data = pd.read_csv(file_path, comment='#', header=None)

        # If merged_data is empty, just copy the first chunk
        if merged_data.empty:
            current_data.columns = columns
            merged_data = current_data.copy()
            cont += 1
        else:
            try:
                # If DataFrame is not named, you can set column names
                current_data.columns = columns

                for col in ['Sum', 'Count_in_Bin', 'Second_Moment', 'Histories_with_Scorer_Active']:
                    if col in current_data.columns:
                        merged_data[col] += current_data[col]

                # Min and Max
                if 'Min' in current_data.columns and 'Min' in merged_data.columns:
                    merged_data['Min'] = np.minimum(merged_data['Min'], current_data['Min'])
                if 'Max' in current_data.columns and 'Max' in merged_data.columns:
                    merged_data['Max'] = np.maximum(merged_data['Max'], current_data['Max'])

                cont += 1
            except:
                pass

    if 'Mean' in columns and 'Sum' in columns and 'Histories_with_Scorer_Active':
        merged_data['Mean'] = merged_data['Sum'] / merged_data['Histories_with_Scorer_Active']
    if 'Mean' in columns and 'Variance' in columns and 'Second_Moment' in columns and 'Histories_with_Scorer_Active':
        merged_data['Variance'] = (merged_data['Second_Moment'] / merged_data['Histories_with_Scorer_Active'] -
                                   merged_data['Mean']**2)
    if 'Standard_Deviation' in columns and 'Variance' in columns:
        merged_data['Standard_Deviation'] = np.sqrt(merged_data['Variance'])


    # Write the combined results to a new output file
    with open(os.path.join(output_path, f'combined_{filename}'), "w") as f:
        for line in header_lines:
            f.write(line)

    merged_data.to_csv(os.path.join(output_path, f'combined_{filename}'), mode='a', header=False, index=False)


    print('')
    print(f'File {filename:}')
    print(f'Number of results merged: {cont} out of {len(output_file_paths)}')
    print('')
