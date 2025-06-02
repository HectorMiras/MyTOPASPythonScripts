import pandas as pd
import re

def read_dnadamage_phase_space(filebase):
    """
    Read the DNADamage phase space data from files with known extensions.
    
    Parameters:
        filebase (str): The base file path/name without extension.
                        The function expects to find filebase+'.header' and filebase+'.phsp'.
    
    Returns:
        pd.DataFrame: The phase space data with proper column names.
    """
    header_file = filebase + ".header"
    phsp_file = filebase + ".phsp"

    # Parse the header file to extract column names.
    with open(header_file, 'r') as f:
        lines = f.readlines()

    columns = []
    header_found = False
    for line in lines:
        line = line.strip()
        if "Columns of data are as follows:" in line:
            header_found = True
            continue
        if header_found:
            if not line:
                continue
            match = re.match(r'\d+:\s*(.+)', line)
            if match:
                colname = match.group(1).strip()
                columns.append(colname)

    # Read the phase-space file using whitespace delimiter.
    df = pd.read_csv(phsp_file, sep='\s+', header=None)
    
    # If there are more column names than columns in the df, slice the header list.
    df.columns = columns[:df.shape[1]]
    
    return df


def total_damage(df):
    # Calculate the total number of particles.
    total_events = len(df)
    
    # Calculate the total energy deposited.
    # Assuming the header file provides energy in keV, convert to MeV.
    total_energy = df['Energy_imparted_per_event [keV]'].sum() / 1000

    # List of damage columns to sum, adjust list as needed.
    damage_columns = [
        'DSBs', 'DSBs_Direct', 'DSBs_Indirect', 'DSBs_Hybrid',
        'DSBs_Direct_WithOneQuasiDirect', 'DSBs_Direct_WithBothQuasiDirect', 'DSBs_Hybrid_WithOneQuasiDirect',
        'SSBs', 'SSBs_Direct', 'SSBs_QuasiDirect', 'SSBs_Indirect',
        'SBs', 'SBs_Direct', 'SBs_QuasiDirect', 'SBs_Indirect',
        'SSB+s', 'DSB+s', 'More complex damages', 'BDs', 'BDs_Direct', 'BDs_QuasiDirect', 'BDs_Indirect',
        'Foci_150nm', 'Foci_500nm'
    ]
    
    damage_totals = {}
    for col in damage_columns:
        if col in df.columns:  # ensure the column exists in the dataframe
            damage_totals[col] = df[col].sum()
    
    # Print the results.
    print(f"Total number of events: {total_events}")
    print(f"Total energy deposited: {total_energy:.2f} MeV")
    for damage, total in damage_totals.items():
        print(f"Total {damage}: {total}")


def merge_dnadamage_files(filebases):
    """
    Merge DNADamage phase space DataFrames from multiple simulation runs.
    
    Parameters:
        filebases (list[str]): List of base file paths/names without extension.
                               Each base file is expected to have corresponding '.header' and '.phsp' files.
    
    Returns:
        pd.DataFrame: A merged DataFrame containing data from all simulation runs.
    """
    dataframes = []
    for base in filebases:
        df = read_dnadamage_phase_space(base)
        dataframes.append(df)
    
    merged_df = pd.concat(dataframes, ignore_index=True)
    return merged_df

