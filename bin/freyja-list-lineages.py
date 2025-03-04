#!/usr/bin/env python3

# Import modules
import pandas as pd
import sys

# written by Adeliza Realingo (03 March 2025)

# Input sample name, tsv file from freyja_demix
sample_tsv = sys.argv[1]
freyja_demix = sys.argv[2]

# Open freyja demix tsv as pandas dataframe
df_sample_lineages = pd.read_csv( freyja_demix, sep='\t')

# Rename the first column into 'info'
df_sample_lineages.rename(columns={df_sample_lineages.columns[0]: "info"}, inplace=True)

# Extract the lineages row
lineages_row = df_sample_lineages.loc[df_sample_lineages['info'] == 'lineages'].values[0]

# Create a new df with the lineages
df_lineages = pd.DataFrame(lineages_row, columns=['lineage'])
df_lineages = df_lineages.drop([0]) # Drops the firs row containing 'lineages'
df_lineages = df_lineages['lineage'].str.split(' ').explode().reset_index(drop=True) # Splits the values by the space inbetween them, then explodes to create new rows for each item


# Write the df to file
df_lineages.to_csv( sample_tsv, sep='\t', index=False )
