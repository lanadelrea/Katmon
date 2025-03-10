#!/usr/bin/env python3

import pandas as pd
import sys

pangolin_csv = sys.argv[1]
nextclade_tsv = sys.argv[2]
sample = sys.argv[3]

# Read the pangolin csv and nextclade tsv files
pangolin = pd.read_csv(pangolin_csv, sep=',')
nextclade = pd.read_csv(nextclade_tsv, sep='\t')


# Clean the dataframes
pangolin.rename(columns={'taxon':'Sample'}, inplace=True) # Rename 'taxon' column as 'sample' in pangolin csv
nextclade.drop(columns=['index'], inplace=True) # Drop 'index' column from nextclade tsv
nextclade.rename(columns={'seqName':'Sample'}, inplace=True) # Rename 'seqName' column as 'sample'

# Merge the dataframes
merged_df = pd.merge(pangolin, nextclade, on='Sample', how='inner')

# Keep only neccessary columns
columns_to_keep = ['Sample', 'lineage', 'scorpio_call', 'note', 'Nextclade_pango', 'clade_who', 'coverage']
lineage_assignment_df = merged_df.loc[:, columns_to_keep].copy()

# Rearrange and rename columns
note_column = lineage_assignment_df.pop('note')
lineage_assignment_df['note'] = note_column
coverage_column = lineage_assignment_df.pop('coverage')
lineage_assignment_df.insert(1, 'coverage', coverage_column)
lineage_assignment_df.rename(columns={'coverage':'Coverage','lineage':'Pangolin', 'scorpio_call': 'Scorpio', 'Nextclade_pango': 'Nextclade', 'clade_who': 'WHO clade', 'note':'Note'}, inplace=True)
lineage_assignment_df.iloc[:, 1] = (lineage_assignment_df.iloc[:, 1] * 100).round(2).astype(str) + '%' 

# Write the dataframe into a new tsv file
lineage_assignment_df.to_csv(f'{sample}_lineage_assignment_merged.tsv', sep='\t', index=False)