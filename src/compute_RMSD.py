from helper_functions import *
import argparse
import numpy as np
import os

# create the argument parser
parser = argparse.ArgumentParser(description='Description of your script')

# add arguments
parser.add_argument('-s', '--step', dest='step', type=int, help='Step for conformations sample extraction')
parser.add_argument('-f', '--file', dest='filename', help='Output file name')

# parse arguments
args = parser.parse_args()

step = args.step
filename = args.filename

# get the current working directory
cwd = os.getcwd()

# navigate to the data directory
data_dir = os.path.join(cwd, 'data', 'aladip')
os.chdir(data_dir)

# Read the entire aladip implicit
print(">>> Reading the raw file ...")
aladip_implicit = read_xyz("aladip_implicit.xyz", 14207380, 3)

# Get conformations and get subsample
print(">>> Processing raw data and extracting a sample of conformations ...")
conformations = get_conformations(aladip_implicit)
conform_sample = get_sample_conformation(conformations, step)

# Compute the RMSD matrix for this sample conformations
print(">>> Computing the RMSD matrix for the sample of conformations...")
print(len(conform_sample))
RMSD_m_sample = compute_RMSD_matrix(conform_sample)

# navigate to the data directory
output_dir = os.path.join(cwd, 'data')
os.chdir(output_dir)

## Save the matrix in a text file
print(">>> Saving the RMSD matrix ...")
np.savetxt(f"{filename}", RMSD_m_sample)

print("Job Done Sucessfully !")