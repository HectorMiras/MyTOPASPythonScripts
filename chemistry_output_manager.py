import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt

def read_irtgvalue_phase_space(filebase):
    """
    Read IRTGValue phase space data from file.
    
    Parameters:
        filebase (str or pathlib.Path): The base file path/name without extension.
                        The function expects filebase+'.phsp' to exist.
    
    Returns:
        pd.DataFrame: The phase space data with proper column names, or None if file doesn't exist.
    """
    # Convert to Path object if it's not already
    filebase = pathlib.Path(filebase)
    
    # Use with_suffix() or with_name() methods for Path objects
    phsp_file = filebase.with_suffix('.phsp')
    
    # Check if file exists
    if not phsp_file.exists():
        print(f"Warning: File {phsp_file} does not exist")
        return None
    
    # Column definitions (per header)
    # 1: GValue       (molecules / 100 eV)
    # 2: GValue_err   (1 σ)
    # 3: Time_ps      (picoseconds)
    # 4: Molecule     (string)
    cols = ['GValue', 'GValue_err', 'Time_ps', 'Molecule']
    
    # Read the phase-space file
    df = pd.read_csv(
        phsp_file,
        comment='#',         # ignore lines starting with '#'
        sep=r'\s+',          # whitespace-delimited
        names=cols,
        header=None,
        engine='python'
    )
    
    return df

def filter_inactive_species(df):
    """
    Filter out species whose total production is zero.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing IRTGValue data.
    
    Returns:
        pd.DataFrame: Filtered DataFrame with only active species.
    """
    active_species = (
        df.groupby('Molecule', sort=False)['GValue']
          .sum()                 # total across all times and histories
          .loc[lambda s: s > 0]  # keep only > 0
          .index
    )
    
    return df[df['Molecule'].isin(active_species)].copy()

def create_pivots(df):
    """
    Create pivot tables for GValue and GValue_err.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing IRTGValue data.
    
    Returns:
        tuple: (pivot, pivot_err) DataFrames with Time_ps as index and Molecule as columns.
    """
    pivot = df.pivot_table(
        index='Time_ps', 
        columns='Molecule',
        values='GValue', 
        aggfunc='mean'
    )
    
    pivot_err = df.pivot_table(
        index='Time_ps', 
        columns='Molecule',
        values='GValue_err', 
        aggfunc='mean'
    )
    
    return pivot, pivot_err

def get_final_values(pivot, pivot_err):
    """
    Extract values at the final time point.
    
    Parameters:
        pivot (pd.DataFrame): Pivot table with Time_ps as index and Molecule as columns.
        pivot_err (pd.DataFrame): Pivot table of errors with same structure.
    
    Returns:
        tuple: (summary, summary_err) Series with final values and errors.
    """
    last_time = pivot.index.max()
    summary = pivot.loc[last_time].to_frame(name='GValue_final')
    summary_err = pivot_err.loc[last_time]
    
    return summary, summary_err, last_time

def plot_final_values(summary, summary_err, last_time=None):
    """
    Create a horizontal bar chart of final values with error bars.
    
    Parameters:
        summary (pd.DataFrame): DataFrame with final GValues.
        summary_err (pd.Series): Series with error values.
        last_time (float, optional): The time of the final values for the title.
    """
    plt.figure(figsize=(6, 4))
    plt.barh(
        summary.sort_values('GValue_final').index,
        summary.sort_values('GValue_final')['GValue_final'].values,
        xerr=summary_err.loc[summary.sort_values('GValue_final').index].values,
        align='center',
        ecolor='k',
        capsize=3
    )
    
    title = 'Final chemical production per species'
    if last_time is not None:
        plt.xlabel(f'GValue (molecules / 100 eV) at t = {last_time:.1f} ps')
    else:
        plt.xlabel('GValue (molecules / 100 eV)')
    
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_time_evolution(pivot, pivot_err):
    """
    Create a time evolution plot with uncertainty bands.
    
    Parameters:
        pivot (pd.DataFrame): Pivot table with Time_ps as index and Molecule as columns.
        pivot_err (pd.DataFrame): Pivot table of errors with same structure.
    """
    plt.figure(figsize=(8, 5))
    for mol in pivot.columns:
        plt.plot(pivot.index, pivot[mol], label=mol)
        plt.fill_between(
            pivot.index,
            pivot[mol] - pivot_err[mol],
            pivot[mol] + pivot_err[mol],
            alpha=0.2
        )
    
    plt.xscale('log')
    plt.xlabel('Time (ps)')
    plt.ylabel('GValue (molecules / 100 eV)')
    plt.title('Temporal evolution of chemical species')
    plt.legend(fontsize='small')
    plt.tight_layout()
    plt.show()

def save_results(summary, pivot, output_dir):
    """
    Save summary and time evolution data to CSV files.
    
    Parameters:
        summary (pd.DataFrame): DataFrame with final GValues.
        pivot (pd.DataFrame): Pivot table with time evolution data.
        output_dir (str or pathlib.Path): Directory to save files.
    """
    output_dir = pathlib.Path(output_dir)
    summary.to_csv(output_dir / 'IRTGValue_GValueFinal.csv')
    pivot.to_csv(output_dir / 'IRTGValue_GValue_vs_time.csv')
    print('✔  Summary and time-evolution CSVs saved.')

def calculate_enhancement(base_filebase, np_filebase):
    """
    Calculate enhancement factors between two simulations (e.g., with/without nanoparticles).
    
    Parameters:
        base_filebase (str or pathlib.Path): Base filepath for baseline simulation.
        np_filebase (str or pathlib.Path): Base filepath for nanoparticle simulation.
    
    Returns:
        tuple: (ratio, ratio_err, t_final) Enhancement factors, errors, and final time.
    """
    # Load data
    df_base = read_irtgvalue_phase_space(base_filebase)
    df_np = read_irtgvalue_phase_space(np_filebase)
    
    # Create pivot tables
    pv_base, pv_base_err = create_pivots(df_base)
    pv_np, pv_np_err = create_pivots(df_np)
    
    # Find common final time
    t_final = min(pv_base.index.max(), pv_np.index.max())
    
    # Get values at final time
    G_base = pv_base.loc[t_final]
    G_np = pv_np.loc[t_final]
    σ_base = pv_base_err.loc[t_final]
    σ_np = pv_np_err.loc[t_final]
    
    # Calculate enhancement factor and error
    ratio = G_np / G_base.replace(0, np.nan)  # avoid /0 warnings
    ratio_err = ratio * np.sqrt((σ_np/G_np)**2 + (σ_base/G_base)**2)
    
    # Remove species with zero production in either simulation
    valid = (G_np > 0) & (G_base > 0)
    ratio = ratio[valid]
    ratio_err = ratio_err[valid]
    
    return ratio, ratio_err, t_final

def plot_enhancement(ratio, ratio_err, t_final=None):
    """
    Create a horizontal bar chart of enhancement factors with error bars.
    
    Parameters:
        ratio (pd.Series): Series with enhancement factors.
        ratio_err (pd.Series): Series with error values.
        t_final (float, optional): The time of enhancement for the title.
    """
    ratio_sorted = ratio.sort_values()
    err_sorted = ratio_err[ratio_sorted.index]
    
    plt.figure(figsize=(7, 4))
    plt.barh(
        ratio_sorted.index, 
        ratio_sorted.values,
        xerr=err_sorted.values, 
        align='center', 
        ecolor='k', 
        capsize=3
    )
    
    plt.axvline(1, ls='--', lw=1, color='grey')
    plt.xlabel('Chemical-species enhancement factor (AuNP / baseline)')
    
    if t_final is not None:
        plt.title(f'Final enhancement (t = {t_final:.1f} ps)')
    else:
        plt.title('Chemical species enhancement')
    
    plt.tight_layout()
    plt.show()

def save_enhancement(ratio, ratio_err, output_path):
    """
    Save enhancement data to a CSV file.
    
    Parameters:
        ratio (pd.Series): Series with enhancement factors.
        ratio_err (pd.Series): Series with error values.
        output_path (str or pathlib.Path): Path to save the CSV file.
    """
    out = pd.DataFrame({
        'Enhancement': ratio,
        'Enhancement_err': ratio_err
    })
    out.to_csv(output_path)
    print(f'✔  CSV "{output_path}" written.')

def merge_multi_run_irtgvalue(base_dir, run_template, n_runs, phsp_name='IRTGValue.phsp'):
    """
    Aggregate data from multiple runs of the same simulation.
    
    Parameters:
        base_dir (str or pathlib.Path): Base directory containing run folders.
        run_template (str): Template string for run folder names, e.g., 'run{}'
        n_runs (int): Number of runs to aggregate.
        phsp_name (str, optional): Name of phase space file. Default is 'IRTGValue.phsp'.
    
    Returns:
        tuple: (df_all, pivot_mean, pivot_sd, summary_df) Various aggregated data.
    """
    base_dir = pathlib.Path(base_dir)
    out_dir = base_dir / 'results'
    out_dir.mkdir(exist_ok=True)
    
    cols = ['GValue', 'GValue_err', 'Time_ps', 'Molecule']
    
    # Load every run into a single MultiIndex DataFrame
    all_runs = []
    for r in range(1, n_runs + 1):
        phsp_path = base_dir / run_template.format(r) / phsp_name
        if not phsp_path.is_file():
            raise FileNotFoundError(f'{phsp_path} not found')
        df = pd.read_csv(
            phsp_path, 
            comment='#', 
            sep=r'\s+',
            names=cols, 
            header=None, 
            engine='python'
        )
        df['Run'] = r
        all_runs.append(df)
    
    df_all = pd.concat(all_runs, ignore_index=True)
    
    # Aggregate time-molecule matrices (mean and SD across runs)
    pivot_mean = (df_all
                  .groupby(['Time_ps', 'Molecule'])['GValue']
                  .mean()
                  .unstack('Molecule')
                  .sort_index())
    
    pivot_sd = (df_all
                .groupby(['Time_ps', 'Molecule'])['GValue']
                .std(ddof=1)
                .unstack('Molecule')
                .reindex_like(pivot_mean))
    
    # Final-time summary (mean ± SD over runs)
    # Collect the final line of each run, then aggregate
    final_rows = []
    for r in range(1, n_runs + 1):
        df_r = df_all[df_all['Run'] == r]
        last_t = df_r['Time_ps'].max()
        final_rows.append(df_r[df_r['Time_ps'] == last_t])
    
    df_final = pd.concat(final_rows)
    
    summary_mean = (df_final
                    .groupby('Molecule')['GValue']
                    .mean()
                    .rename('Mean_GValue'))
    
    summary_sd = (df_final
                  .groupby('Molecule')['GValue']
                  .std(ddof=1)
                  .rename('SD_GValue'))
    
    # Remove zero-production species
    nonzero = summary_mean > 0
    summary_mean = summary_mean[nonzero]
    summary_sd = summary_sd[nonzero]
    
    summary_df = pd.concat([summary_mean, summary_sd], axis=1)
    
    # Save aggregates to results
    pivot_mean.to_csv(out_dir / 'IRTGValue_GValue_vs_time_mean.csv')
    pivot_sd.to_csv(out_dir / 'IRTGValue_GValue_vs_time_sd.csv')
    summary_df.to_csv(out_dir / 'IRTGValue_GValue_final_mean_sd.csv')
    
    # Create per-run table at final time
    pivot_final = (df_final
                   .pivot(index='Run', columns='Molecule', values='GValue')
                   .reindex(columns=summary_mean.index))
    
    pivot_final.to_csv(out_dir / 'IRTGValue_GValue_final_per_run.csv')
    
    print(f'✔  Aggregated files written to {out_dir}')
    
    return df_all, pivot_mean, pivot_sd, summary_df, pivot_final, out_dir

def plot_multi_run_summary(summary_mean, summary_sd):
    """
    Create a bar chart of final values with error bars from multiple runs.
    
    Parameters:
        summary_mean (pd.Series): Series with mean values.
        summary_sd (pd.Series): Series with standard deviation values.
    """
    plt.figure(figsize=(7, 4))
    summary_mean_sorted = summary_mean.sort_values()
    plt.barh(
        summary_mean_sorted.index,
        summary_mean_sorted.values,
        xerr=summary_sd.loc[summary_mean_sorted.index].values,
        align='center', 
        ecolor='k', 
        capsize=3
    )
    
    plt.xlabel('Mean GValue (molecules / 100 eV)')
    plt.title('Final chemical production per species (multiple runs)')
    plt.tight_layout()
    plt.show()

def plot_multi_run_violin(pivot_final):
    """
    Create a violin plot showing distribution of values across runs.
    
    Parameters:
        pivot_final (pd.DataFrame): DataFrame with Run as index and Molecule as columns.
    """
    species = pivot_final.columns
    data = [pivot_final[col].values for col in species]
    
    plt.figure(figsize=(max(6, 0.6*len(species)), 4))
    plt.violinplot(data, showmeans=True, showmedians=False, showextrema=False)
    
    plt.xticks(np.arange(1, len(species)+1), species, rotation=45, ha='right')
    plt.ylabel('Final GValue (molecules / 100 eV)')
    plt.title('Distribution of final chemical production across runs')
    plt.tight_layout()
    plt.show()