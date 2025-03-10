#!/usr/bin/env python3

# written by Adeliza Realingo (07 March 2025)

# Import modules
import pandas as pd
import os
import re
import sys

# Input: Folder path containing two files
folder_path = sys.argv[1]

# List all files in the folder
files = os.listdir(folder_path)

# Ensure there are exactly two files in the folder
if len(files) != 2:
    raise ValueError("The folder must contain exactly two files.")

# Assign each file to A_lineage_mutations_path and B_lineage_mutations_path
A_lineage_mutations_path = os.path.join(folder_path, files[0])
B_lineage_mutations_path = os.path.join(folder_path, files[1])

# Input
#A_lineage_mutations_path = sys.argv[1]
#B_lineage_mutations_path = sys.argv[2]
sample = sys.argv[2]

df_A_lineage = pd.read_csv( A_lineage_mutations_path, sep="\t")
df_B_lineage = pd.read_csv( B_lineage_mutations_path, sep="\t")

# Get value for first (A) and second (B) lineage, then use the value as column name
A_lineage = (os.path.basename(A_lineage_mutations_path).split('_')[-1]).rsplit('.', 1)[0]
A_lineage_mutations = pd.read_csv( A_lineage_mutations_path, sep='\t', header=None)
A_lineage_mutations.columns = [A_lineage]
# ~ #
B_lineage = (os.path.basename(B_lineage_mutations_path).split('_')[-1]).rsplit('.', 1)[0]
B_lineage_mutations = pd.read_csv( B_lineage_mutations_path, sep='\t', header=None)
B_lineage_mutations.columns = [B_lineage]

# Extract only the nucleotide mutations and convert series
A_mutations_df = A_lineage_mutations[A_lineage].str.extract(r'([A-Za-z0-9]+)')
A_mutations_ser = A_mutations_df.iloc[:,0]
# ~ #
B_mutations_df = B_lineage_mutations[B_lineage].str.extract(r'([A-Za-z0-9]+)')
B_mutations_ser = B_mutations_df.iloc[:,0]

# Extract the position and mutation
A_extracted_pos_mut = A_mutations_ser.apply(lambda x: tuple(map(lambda y: int(y) if y.isdigit() else y, re.findall(r'(\d+)([A-Za-z])$', x)[0])))
B_extracted_pos_mut = B_mutations_ser.apply(lambda x: tuple(map(lambda y: int(y) if y.isdigit() else y, re.findall(r'(\d+)([A-Za-z])$', x)[0])))

# Format the output
A_output = '\n'.join(f"{num} {letter}" for num, letter in A_extracted_pos_mut)
B_output = '\n'.join(f"{num} {letter}" for num, letter in B_extracted_pos_mut)

# Save to their corresponding text files
with open(f'{sample}_lineage_A.txt', 'w') as f:
    f.write(A_output)

with open(f'{sample}_lineage_B.txt', 'w') as f:
    f.write(B_output)