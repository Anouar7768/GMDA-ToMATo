# Code description 

There are two python scripts which are located in the **src** folder.

1) `compute_RMSD.py` script is the first script that computes the RMSD matrix for a sample of conformations. To run `compute_RMSD.py` script, be sure you are in the folder **GMDA-ToMATo** (not in GMDA-ToMATo/src), and specify the following arguments:
- **step** : the step at which we bypass the conformations, in order to create a significative sample of conformations, and avoid selecting only neighbors.
- **filename** : the filename of the output file where we will store the computed RMSD matrix. By default, this output is sotred in the data directory. The output is a text file, so you need to specifiy the .txt extension

*Example* : `python src/compute_RMSD.py -s 400 -f test.txt`  

2)`run_ToMATo.py` is the second script that runs the ToMATo algorithm and store the labels of clustering. To run `run_ToMATo.py` script, be sure you are in the folder **GMDA-ToMATo** (not in GMDA-ToMATo/src), and specify the following arguments:
- **in_filename** the filename of the input file, where the RMSD matrix is stored.
- **out_filename** the filename of the output file, where we will store the labels
- **n_clusters** the number of clusters for the ToMATo algorithm 

*Example* : `python src/run_ToMATo.py -if 400 test.txt -of output.txt -n 5`  
