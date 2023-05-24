import struct
import pandas as pd
import numpy as np
import os
from phsp_manager import generate_header


# Read header file
def read_header(header_file):

    iaea_header_params = {}

    with open(header_file, 'r') as file:
        header_lines = file.readlines()

    recordlength = 0

    for line in header_lines:
        if "// X is stored ?" in line:
            iaea_header_params['X'] = int(line.split()[0])
        elif "// Y is stored ?" in line:
            iaea_header_params['Y'] = int(line.split()[0])
        elif "// Z is stored ?" in line:
            iaea_header_params['Z'] = int(line.split()[0])
        elif "// U is stored ?" in line:
            iaea_header_params['U'] = int(line.split()[0])
        elif "// V is stored ?" in line:
            iaea_header_params['V'] = int(line.split()[0])
        elif "// W is stored ?" in line:
            iaea_header_params['W'] = int(line.split()[0])
        elif "// Weight is stored ?" in line:
            iaea_header_params['Weight'] = int(line.split()[0])
        elif "// Extra floats stored ?" in line:
            iaea_header_params['ExtraFloats'] = int(line.split()[0])
        elif "// Extra longs stored ?" in line:
            iaea_header_params['ExtraLongs'] = int(line.split()[0])
        elif "// Constant X" in line:
            iaea_header_params['Xconst'] = float(line.split()[0])
        elif "// Constant Y" in line:
            iaea_header_params['Yconst'] = float(line.split()[0])
        elif "// Constant Z" in line:
            iaea_header_params['Zconst'] = float(line.split()[0])

    recordlength = 4 * (iaea_header_params['X'] + iaea_header_params['Y'] + iaea_header_params['Z'])
    recordlength += 4 * (iaea_header_params['U'] + iaea_header_params['V'] + iaea_header_params['Weight'])
    recordlength += iaea_header_params['W']
    recordlength += 4 * (iaea_header_params['ExtraFloats'] + iaea_header_params['ExtraLongs'])
    recordlength += 4  # Particle Energy
    iaea_header_params['RECORD_LENGTH'] = recordlength

    return iaea_header_params


def read_particle(bytearr, iaea_header_params):
    particle = {}
    Wsign = 1
    IsNewHistory = False
    kpar_dic = {1:22, 2:11, 3:-11}
    # Particle type
    kpar = struct.unpack('b', bytearr[:1])[0]
    if kpar < 0:
        Wsign = -1
        kpar = -kpar

    particle['kpar'] = kpar_dic[kpar]

    # Energy
    E = struct.unpack('f', bytearr[1:5])[0]
    IsNewHistory = E < 0
    particle['E'] = abs(E)

    offset = 5  # bytearr index offset
    parameters = ['X', 'Y', 'Z', 'U', 'V', 'Weight']
    for param in parameters:
        if iaea_header_params[param] == 1:
            particle[param] = struct.unpack('f', bytearr[offset:offset+4])[0]
            offset += 4
        else:
            particle[param] = 0  # or some constant, as in C# code

    # Compute W if needed
    if iaea_header_params['W'] == 1:
        aux = particle['U']**2 + particle['V']**2
        if aux <= 1.0:
            particle['W'] = Wsign * (1 - aux)**0.5
        else:
            aux = aux**0.5
            particle['U'] /= aux
            particle['V'] /= aux

    # Read extra floats
    extra_floats = []
    for i in range(iaea_header_params['ExtraFloats']):
        extra_floats.append(struct.unpack('f', bytearr[offset:offset + 4])[0])
        offset += 4
    for i, value in enumerate(extra_floats):
        particle[f'extra_float{i + 1}'] = value

    # Read extra longs
    extra_longs = []
    for i in range(iaea_header_params['ExtraLongs']):
        extra_longs.append(struct.unpack('i', bytearr[offset:offset + 4])[0])
        offset += 4
    for i, value in enumerate(extra_longs):
        particle[f'extra_long{i + 1}'] = value

            
    particle['IsNewHistory'] = IsNewHistory

    return particle


# Write data to file in TOPAS format
def write_topas_file(data, filename):
    data.to_csv(filename, sep=' ', index=False, header=False)


# Transform IAEA to TOPAS
def transform_iaea_to_topas(header_file, binary_file, topas_file):
    # Read the header file
    iaea_header_params = read_header(header_file)

    # Read the binary file
    with open(binary_file, 'rb') as f:
        byte_data = f.read()

    # Prepare an empty dataframe for the particle data
    columns = ['kpar', 'IsNewHistory', 'E', 'X', 'Y', 'Z', 'U', 'V', 'W', 'Weight']
    data = pd.DataFrame(columns=columns)

    record_length = iaea_header_params['RECORD_LENGTH']
    num_records = len(byte_data) // record_length

    #data_list = []
    with open(topas_file, 'w') as f:
        for i in range(num_records):
            record = byte_data[i * record_length:(i + 1) * record_length]
            particle = read_particle(record, iaea_header_params)
            particle_line = "{:14.9f} {:14.9f} {:14.9f} {:14.9f} {:14.9f} {:14.9f} {:14.9f} {:14d} {:d} {:d}".format(
                particle['X'], particle['Y'], particle['Z'], particle['U'], particle['V'], particle['E'],
                particle['Weight'], particle['kpar'], int(particle['W']), particle['IsNewHistory'])
            f.write(particle_line + "\n")
            #particle_df = pd.DataFrame(particle, index=[0])
            #data_list.append(particle_df)

    # Write the data to the TOPAS file
    #data = pd.concat(data_list, ignore_index=True)
    #write_topas_file(data, topas_file)


# Call the function
folder_path = '../tests/phsp_IAEA2TOPAS/'
header_file = 'TrueBeam6FFF_10x10.IAEAheader'
binary_file = 'TrueBeam6FFF_10x10.IAEAphsp'
topas_file = 'TrueBeam6FFF_10x10'
transform_iaea_to_topas(folder_path+header_file, folder_path+binary_file, folder_path+topas_file+'.phsp')
generate_header(folder_path, topas_file)
