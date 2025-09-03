#!/usr/bin/env python3

import sys
import re
import os
import pandas as pd
import subprocess as sp

sample_name = sys.argv[1]
bammix_csv = sys.argv[2]

# Set threshold
mix_thresh = float(sys.argv[3]) # This is the proportion of the major allele
#mix_thresh = 0.8
depth_thresh = 20 # Read depth
pos_thresh = 4 # Minimum number of positions with mixtures to flag by bammix

df = pd.read_csv(bammix_csv) # Put the csv file into a dataframe

# Get df of nucleotide positions
df_pos = df[['Position']]
# Get positions with >20% mixture from df
df_prop = df[["A_proportion", "C_proportion", "G_proportion", "T_proportion"]]
# Get largest number for each row in df
df_prop_max = df_prop.max(axis=1)
# Select rows with less than 0.80 value in df_prop_max
df_prop_max_bool = df_prop_max < mix_thresh
df_prop_below80max = df_prop[df_prop_max_bool]
# Flag only if depth is >20 reads
df_depth = df[["Total_reads"]]
df_depth_20x = df_depth["Total_reads"] > depth_thresh
df_prop_below80max_20xdepth = df_prop_below80max[df_depth_20x]
second_max_value = df_prop_below80max_20xdepth.apply(lambda row: row.nlargest(2).iloc[-1], axis=1)
# Flag only if depth is >20 reads, rows have <0.80 value in df_prop_max, and second max value > 0.20
df_flagged = df_prop_below80max_20xdepth[(df_prop_max_bool) & df_depth_20x & (second_max_value > 0.20)]

bad_mixtures_positions = df_pos.loc[df_flagged.index, "Position"].tolist()
bad_mixtures_base_prop = df_flagged.values.tolist()

# Write table of flagged barcodes, positions, base proportions
bammix_bad_mixtures_csv = pd.DataFrame()
bammix_bad_mixtures_csv['barcode'] = sample_name
bammix_bad_mixtures_csv['positions'] = bad_mixtures_positions
bammix_bad_mixtures_csv['base_prop'] = bad_mixtures_base_prop
bammix_bad_mixtures_csv.to_csv(f'{sample_name}_bammix_flags.csv', index=False)