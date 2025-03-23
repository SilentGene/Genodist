import pandas as pd
import numpy as np
import seaborn as sns
import argparse

# usage: python draw_heatmap.py -i ANI_results.tsv -o ani_clustermap.pdf
# This script will ignore the first row of the input file, and only consider the first three columns

"""exzample of input tsv file
Genome1	Genome2	ANI(or AAI)
genome1	genome2	99.9
genome1	genome3	98.7
genome2	genome3	97.8
"""


parser = argparse.ArgumentParser(description='Draw a clustermap of ANI/AAI values')
parser.add_argument('-i', required=True, help='Input ANI/AAI results file. TSV format with 3 columns: First genome, Second genome, ANI/AAI value')
parser.add_argument('-o', required=True, help='Output PDF file')
args = parser.parse_args()

# Read data
data = pd.read_csv(args.i, sep='\t', header=0)

# Get unique genome names from the first two columns (index 0 and 1)
genomes = sorted(set(data.iloc[:, 0]).union(set(data.iloc[:, 1])))

# Create ANI matrix
ani_matrix = pd.DataFrame(index=genomes, columns=genomes, dtype=float)

# Fill ANI matrix using column indices
for _, row in data.iterrows():
    ani_matrix.loc[row.iloc[0], row.iloc[1]] = float(row.iloc[2])

# Fill diagonal with 100%
np.fill_diagonal(ani_matrix.values, 100.0)

# Ensure all values are floats
ani_matrix = ani_matrix.astype(float)

# Create a seaborn clustermap with annotations, and save it as a PDF
clustermap = sns.clustermap(ani_matrix, annot=True, cmap="coolwarm", fmt=".2f", figsize=(12, 12))

# Save the clustermap to a PDF file
clustermap.savefig(args.o)
