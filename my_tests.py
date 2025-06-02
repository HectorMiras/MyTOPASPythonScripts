import dnadamage_phsp_manager 
import sddparser
import pprint


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

# To run the test:
test_parseSDDFile()