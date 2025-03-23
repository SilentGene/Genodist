import pandas as pd
import numpy as np
import seaborn as sns
import argparse

# usage: python draw_heatmap.py -i all_ani_results.tsv -o ani_clustermap.pdf

# ani and aai visualization are not the same, aai is calculated two-way, so will have to fill the matrix symmetrically

parser = argparse.ArgumentParser(description='Draw a clustermap of Average Nucleotide Identity (AAI) values')
parser.add_argument('-i', required=True, help='Input AAI results file. TSV format with columns: Genome1, Genome2, AAI')
parser.add_argument('-o', required=True, help='Output PDF file')
args = parser.parse_args()

# Read data
data = pd.read_csv(args.i, sep='\t', header=0)
# Get unique genome names
genomes = sorted(set(data['Genome1']).union(set(data['Genome2'])))

# Create AAI matrix
aai_matrix = pd.DataFrame(index=genomes, columns=genomes, dtype=float)

# Fill ANI matrix
for _, row in data.iterrows():
    aai_matrix.loc[row['Genome1'], row['Genome2']] = float(row['AAI'])
    aai_matrix.loc[row['Genome2'], row['Genome1']] = float(row['AAI'])  # Symmetric fill

# Fill diagonal with 100%
np.fill_diagonal(aai_matrix.values, 100.0)

# Ensure all values are floats
aai_matrix = aai_matrix.astype(float)

# Create a seaborn clustermap with annotations, and save it as a PDF
clustermap = sns.clustermap(aai_matrix, annot=True, cmap="coolwarm", fmt=".2f", figsize=(12, 12))

# Save the clustermap to a PDF file
clustermap.savefig(args.o)
