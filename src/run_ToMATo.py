import numpy as np
import argparse
import os

from gudhi.clustering.tomato import Tomato


# Setting hyper-parameters
density = "logDTM"
k = 30

# Reading the RMSD matrix
parser = argparse.ArgumentParser(description='Description of your script')

# add arguments
parser.add_argument('-if', '--in_file', dest='in_filename', help='Input file name (.txt)')
parser.add_argument('-of', '--out_file', dest='out_filename', help='Output file name (.txt)')
parser.add_argument('-n', '--n_clusters', dest='n_clusters', type=int, help='Number of clusters')

# parse arguments
args = parser.parse_args()

# get the current working directory
cwd = os.getcwd()

# navigate to the data directory
data_dir = os.path.join(cwd, 'data')
os.chdir(data_dir)

# Read the matrix from the text file
print(">>> Reading the RMSD matrix...")
RMSD_matrix = np.loadtxt(args.in_filename)

print(">>> Running ToMATo Clustering...")
# Create Tomato object
t = Tomato(density_type=density, k=k, metric="precomputed", n_jobs=-1)

# Fit the data
t.fit(RMSD_matrix)
t.plot_diagram()

# Create clusters
t.n_clusters_ = args.n_clusters
labels = t.labels_

print(">>> Saving the labels...")
np.savetxt(f"{args.out_filename}", labels)
