from helper_functions import *
import argparse
import numpy as np
import os

# create the argument parser
parser = argparse.ArgumentParser(description='Description of your script')

# add arguments
parser.add_argument('-k1', '--lower', dest='k1', type=int, help='Lower Bound for the conformation sample extraction')
parser.add_argument('-k2', '--upper', dest='k2', type=int, help='Upper Bound for the conformation sample extraction')
parser.add_argument('-f', '--file', dest='filename', help='Output file name')

# parse arguments
args = parser.parse_args()

k1 = args.k1
k2 = args.k2
filename = args.filename

# get the current working directory
cwd = os.getcwd()

# navigate to the data directory
data_dir = os.path.join(cwd, 'data', 'aladip')
os.chdir(data_dir)

# Read the entire aladip implicit
print(">>> Reading the raw file ...")
aladip_implicit = read_xyz("aladip_implicit.xyz", k2*10, 3)

# Get conformations and get subsample
print(">>> Processing raw data and extracting a sample of conformations ...")
conformations = get_conformations(aladip_implicit)
conform_sample = get_sample_conformation(conformations, k1, k2)

# Compute the RMSD matrix for this sample conformations
print(">>> Computing the RMSD matrix for the sample of conformations...")
RMSD_m_sample = compute_RMSD_matrix(conform_sample)

# navigate to the data directory
output_dir = os.path.join(cwd, 'data')
os.chdir(output_dir)

## Save the matrix in a text file
print(">>> Saving the RMSD matrix ...")
np.savetxt(f"{filename}", RMSD_m_sample)

print("Job Done Sucessfully !")