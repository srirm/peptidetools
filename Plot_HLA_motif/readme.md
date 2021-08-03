# Plot HLA motif from 9mer cores
##Will generate and plot position specific frequencies of amino acids in 9mer peptides

- This will take text file containing list of peptides in one column. Specify this file at start of script

- separates out the aminoacids into position specific columns (Pos- 1-9) into a dataframe

- calculates the number of amino acids at each position and its percentage.
- plots a heatmap

Example file DR_cores.txt has 9-mer minimal cores. 
Output will be stores in folder 'output'

Requires: pandas, seaborn, matplotlib