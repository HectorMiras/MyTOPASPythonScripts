import sys
import os.path
import glob

searchpath = "/home/radiofisica/hector/mytopassimulations/TOPAS_CellsNPs/work/CellColony-med5-cell5"  # Set your base directory here
phspfilename = "PhaseSpace_NP"  # Set your phsp file base name here (without extension)


def search_run_phsps(phsp_name, cell_path):
    """Find all run phsp files in a cell directory."""
    run_dirs = sorted(glob.glob(os.path.join(cell_path, "run*")))
    phsp_files = []
    for run_dir in run_dirs:
        phsp_file = os.path.join(run_dir, f"{phsp_name}.phsp")
        header_file = os.path.join(run_dir, f"{phsp_name}.header")
        if os.path.exists(phsp_file) and os.path.exists(header_file):
            phsp_files.append((phsp_file, header_file))
    return phsp_files


def search_cell_phsps(phsp_name, searchpath):
    """Find all combined run phsp files in all cell directories."""
    cell_dirs = sorted(glob.glob(os.path.join(searchpath, "cell*")))
    phsp_files = []
    for cell_dir in cell_dirs:
        combined_file = os.path.join(cell_dir, "results", f"{phsp_name}_combined_runs.phsp")
        header_file = os.path.join(cell_dir, "results", f"{phsp_name}_combined_runs.header")
        if os.path.exists(combined_file) and os.path.exists(header_file):
            phsp_files.append((combined_file, header_file))
    return phsp_files


def merge_phsp(phsp_list, output_path):
    """Merge phsp files and headers in phsp_list, write to output_path (tuple: (phsp, header)).
    Automatically detects units and the presence of Subcomponent Index column from the first header file."""
    # Initialize accumulators
    hist, hist_reac, scored, ne, ng = 0, 0, 0, 0, 0
    mine, ming, maxe, maxg = [], [], [], []
    noempty = 0
    # Defaults in case header is missing
    pos_unit = "[cm]"
    energy_unit = "[MeV]"
    has_subcomponent_index = False
    # Detect units and subcomponent index from first header
    for _, header_file in phsp_list:
        try:
            with open(header_file, "r") as f:
                for line in f:
                    if line.strip().startswith("1: Position X"):
                        if "[" in line and "]" in line:
                            pos_unit = line[line.find("["):line.find("]")+1]
                    if line.strip().startswith("6: Energy"):
                        if "[" in line and "]" in line:
                            energy_unit = line[line.find("["):line.find("]")+1]
                    if line.strip().startswith("11: Subcomponent Index"):
                        has_subcomponent_index = True
                break
        except Exception:
            continue
    with open(output_path[0], "wb") as out_phsp:
        for phsp_file, header_file in phsp_list:
            try:
                with open(header_file, "r") as f:
                    mensaje = f.readlines()
                    if len(mensaje) > 0:
                        # Helper to extract value from a line starting with a keyword
                        def extract_value(lines, keyword, typ=int, default=0):
                            for line in lines:
                                if line.strip().startswith(keyword):
                                    try:
                                        return typ(line.split(":")[-1].strip().split()[0])
                                    except Exception:
                                        return default
                            return default
                        hist += extract_value(mensaje, "Number of Original Histories", int, 0)
                        hist_reac += extract_value(mensaje, "Number of Original Histories that Reached Phase Space", int, 0)
                        scored += extract_value(mensaje, "Number of Scored Particles", int, 0)
                        ne += extract_value(mensaje, "Number of e-", int, 0)
                        ng += extract_value(mensaje, "Number of gamma", int, 0)
                        # For min/max energies, use float and allow for missing
                        def extract_float(lines, keyword):
                            for line in lines:
                                if line.strip().startswith(keyword):
                                    try:
                                        return float(line.split(":")[-1].strip().split()[0])
                                    except Exception:
                                        return None
                            return None
                        min_e = extract_float(mensaje, "Minimum Kinetic Energy of e-")
                        if min_e is not None:
                            mine.append(min_e)
                        min_g = extract_float(mensaje, "Minimum Kinetic Energy of gamma")
                        if min_g is not None:
                            ming.append(min_g)
                        max_e = extract_float(mensaje, "Maximum Kinetic Energy of e-")
                        if max_e is not None:
                            maxe.append(max_e)
                        max_g = extract_float(mensaje, "Maximum Kinetic Energy of gamma")
                        if max_g is not None:
                            maxg.append(max_g)
                        noempty += 1
                with open(phsp_file, "rb") as f2:
                    data = f2.read()
                    out_phsp.write(data)
            except Exception as e:
                print(f"Problems processing file {phsp_file}: {e}")
    # Write merged header
    with open(output_path[1], "w") as f:
        f.write("TOPAS ASCII Phase Space\n\n")
        f.write(f"Number of Original Histories: {hist}\n")
        f.write(f"Number of Original Histories that Reached Phase Space: {hist_reac}\n")
        f.write(f"Number of Scored Particles: {scored}\n\n")
        f.write("Columns of data are as follows:\n")
        f.write(f" 1: Position X {pos_unit}\n")
        f.write(f" 2: Position Y {pos_unit}\n")
        f.write(f" 3: Position Z {pos_unit}\n")
        f.write(" 4: Direction Cosine X\n")
        f.write(" 5: Direction Cosine Y\n")
        f.write(f" 6: Energy {energy_unit}\n")
        f.write(" 7: Weight\n")
        f.write(" 8: Particle Type (in PDG Format)\n")
        f.write(" 9: Flag to tell if Third Direction Cosine is Negative (1 means true)\n")
        f.write("10: Flag to tell if this is the First Scored Particle from this History (1 means true)\n")
        if has_subcomponent_index:
            f.write("11: Subcomponent Index\n")
        f.write("\n")
        f.write(f"Number of e-: {ne}\n")
        f.write(f"Number of gamma: {ng}\n\n")
        f.write(f"Minimum Kinetic Energy of e-: {min(mine) if mine else 0} {energy_unit}\n")
        f.write(f"Minimum Kinetic Energy of gamma: {min(ming) if ming else 0} {energy_unit}\n\n")
        f.write(f"Maximum Kinetic Energy of e-: {max(maxe) if maxe else 0} {energy_unit}\n")
        f.write(f"Maximum Kinetic Energy of gamma: {max(maxg) if maxg else 0} {energy_unit}\n")


def main():
    # Set this to True if your phsp files have the 11th column
    has_subcomponent_index = False  # Change to True if needed
    # Step 1: Merge runs within each cell
    cell_dirs = sorted(glob.glob(os.path.join(searchpath, "cell*")))
    for cell_dir in cell_dirs:
        phsp_list = search_run_phsps(phspfilename, cell_dir)
        if not phsp_list:
            print(f"No phsp files found in {cell_dir}")
            continue
        results_dir = os.path.join(cell_dir, "results")
        os.makedirs(results_dir, exist_ok=True)
        output_path = (
            os.path.join(results_dir, f"{phspfilename}_combined_runs.phsp"),
            os.path.join(results_dir, f"{phspfilename}_combined_runs.header")
        )
        print(f"Merging {len(phsp_list)} runs in {cell_dir}")
        merge_phsp(phsp_list, output_path)

    # Step 2: Merge all cells
    cell_phsp_list = search_cell_phsps(phspfilename, searchpath)
    if not cell_phsp_list:
        print("No combined cell phsp files found.")
        return
    results_dir = os.path.join(searchpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    output_path = (
        os.path.join(results_dir, f"{phspfilename}_all_cells.phsp"),
        os.path.join(results_dir, f"{phspfilename}_all_cells.header")
    )
    print(f"Merging {len(cell_phsp_list)} cells into {output_path[0]}")
    merge_phsp(cell_phsp_list, output_path)

if __name__ == "__main__":
    main()


