#!/usr/bin/env python3

import argparse

from ast import literal_eval
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pysam
import sys
######################## S E C O N D  P L O T ############################

# Input files
primer_scheme_file = sys.argv[1]                               # SARS-CoV-2.scheme.bed
primer_scheme = pd.read_csv(primer_scheme_file, sep='\t')

read_depth_path = sys.argv[2]
alt_plot = pd.read_csv(read_depth_path, sep='\t')

sample = sys.argv[3]

# Output files
aaf_plot = sys.argv[4]                                         # AAFplot_amplicons.png

############################################################################

# Initialize a dictionary to store primer pair information
primer_pairs = {}

# Iterate through the primer scheme dataframe and fill the primer_pairs dictionary
for index, row in primer_scheme.iterrows():
    primer_id = row['Amplicon_ID']
    primer_number, direction = primer_id.rsplit('_', 1)

    if direction == 'alt1':
        primer_number_without_alt1, alt1_direction = primer_number.rsplit('_', 1)
        if primer_number_without_alt1 not in primer_pairs:
            primer_pairs[primer_number_without_alt1] = {'LEFT': None, 'RIGHT': None, 'ALT1': None}
        if alt1_direction == 'LEFT':
            primer_pairs[primer_number_without_alt1]['ALT1'] = row['Start']
        elif alt1_direction == 'RIGHT':
            primer_pairs[primer_number_without_alt1]['ALT1'] = row['End']
    else:
        if primer_number not in primer_pairs:
            primer_pairs[primer_number] = {'LEFT': None, 'RIGHT': None, 'ALT1': None}
        if direction == 'LEFT':
            primer_pairs[primer_number]['LEFT'] = row['End']
        elif direction == 'RIGHT':
            primer_pairs[primer_number]['RIGHT'] = row['Start']

# Create a list to store the amplicon range data
amplicon_ranges = []

# Calculate and append the amplicon ranges to the list
for primer_number, positions in primer_pairs.items():
    left_end = positions.get('LEFT')
    right_start = positions.get('RIGHT')
    alt1_start = positions.get('ALT1')

    if left_end is not None and right_start is not None:
        amplicon_range = (left_end + 1, right_start - 1)
        amplicon_ranges.append({'Amplicon': primer_number, 'Start': amplicon_range[0], 'End': amplicon_range[1]})
    elif alt1_start is not None:
        amplicon_ranges.append({'Amplicon': primer_number, 'Start': alt1_start, 'End': None})
    else:
        amplicon_ranges.append({'Amplicon': primer_number, 'Start': None, 'End': None})

# Create a DataFrame from the amplicon_ranges list
amplicon_ranges_df = pd.DataFrame(amplicon_ranges)



# Adding annotation
# alt_plot['aa'] = alt_plot.apply(lambda x: f'{x.gene}:{x.protein_change}', axis=1)
# alt_plot['lineage'] = alt_plot.apply(lambda x: get_lineage(x.Delta, x.Omicron), axis=1)
from ast import literal_eval

alt_plot_variant = alt_plot['variant']
alt_plot['variant'] = alt_plot_variant.apply(literal_eval)
alt_plot['variant'] = alt_plot['variant'].apply(lambda x: (int(x[0]), x[1], x[2]))
alt_plot['position'] = alt_plot['variant'].apply(lambda x: int(x[0]))
alt_plot['ref'] = alt_plot['variant'].apply(lambda x: x[1])
alt_plot['alt'] = alt_plot['variant'].apply(lambda x: x[2])

alt_plot['gene'] = alt_plot['aa'].apply(lambda x: x.split(':')[0])
alt_plot['protein_change'] = alt_plot['aa'].apply(lambda x: x.split(':')[1])

### Matching positions dataframe ###

# Create an empty list to store matching positions
matching_positions = []

# Iterate through the mapping DataFrame
for index, row in alt_plot.iterrows():
    position = row['position']
    for amplicon, amplicon_range_row in amplicon_ranges_df.iterrows():
        start = amplicon_range_row['Start']
        end = amplicon_range_row['End']
        if start <= position <= end:
            matching_positions.append({
                'Gene': row['gene'],
                'Position': position,
                'Amplicon': amplicon + 1,
                'Start': start,
                'End': end,
                'Lineage': row['mutation_label'],
                'Protein_change': row['protein_change'],
                'Ref_allele': row['ref'],
                'Alt_allele': row['alt'],
                'Alt_allele_depth': row['alt_ad'],
                'Read_depth': row['dp'],
                'AAF': row['aaf'],
                'Amplicon_num': amplicon + 1
            })

# Create a DataFrame from the list of matching positions
matching_positions_df = pd.DataFrame(matching_positions)

# Set MultiIndex with "Gene" and "Position"
matching_positions_df.set_index(['Gene', 'Amplicon', 'Position'], inplace=True)
#matching_positions_df = matching_positions_df[matching_positions_df['Lineage'].notna()]

mp_df = matching_positions_df

# Group by Amplicon and perform the necessary calculations
grouped = mp_df.groupby(['Amplicon'])

def process_group(group):
    if len(group) == 1:
        row = group.iloc[0]
        if row['Lineage'] == 'Omicron':
            group['Omicron'] = row['AAF']
            group['Delta'] = 1 - row['AAF']
        elif row['Lineage'] == 'Delta':
            group['Delta'] = row['AAF']
            group['Omicron'] = 1 - row['AAF']
        else:
            group['Both'] = 1
    else:
        delta_rows = group[group['Lineage'] == 'Delta']
        omicron_rows = group[group['Lineage'] == 'Omicron']

        if len(delta_rows) == len(group):
            max_aaf = group['AAF'].max()
            group['Delta'] = max_aaf
            group['Omicron'] = 1 - max_aaf
        elif len(omicron_rows) == len(group):
            max_aaf = group['AAF'].max()
            group['Omicron'] = max_aaf
            group['Delta'] = 1 - max_aaf
        else:
            delta_max_aaf = delta_rows['AAF'].max()
            omicron_max_aaf = omicron_rows['AAF'].max()
            group['Delta'] = delta_max_aaf
            group['Omicron'] = omicron_max_aaf

    return group.head(1)  # Keep only the first row with the calculated AAF values

result_df = grouped.apply(process_group).reset_index(drop=True)

#result_df.set_index('Amplicon_num', inplace=True)

# result_df = result_df[result_df['Alt_allele_depth'].notna()]
result_df = result_df[result_df['Alt_allele_depth'] !=0]
result_df = result_df[result_df['Both']!=1]

# Create dataframe with Amplicon, Delta AAF, and Omicron AAF
# df_aaf = df_aaf[df_aaf['Lineage'].notna()]
df_aaf = result_df[["Amplicon_num", "Omicron", "Delta"]].copy()
df_aaf.rename(columns = {'Amplicon_num':'Amplicon'}, inplace=True)
df_aaf.set_index("Amplicon", inplace=True)
df_aaf['Delta'] = 1 - df_aaf['Omicron']

# Create stacked bar chart
df_aaf.plot(kind='bar', stacked=True, color=['hotpink', 'blue'], figsize=(20, 5))

# Labels for x and y axis
plt.xlabel('Amplicon')
plt.ylabel('Alternative Allele Fraction')

# Other plot stuff
plt.title(sample)

plt.xticks(rotation=90)
plt.ylim([0,1])
plt.legend(bbox_to_anchor=[1,1], ncol=1)
plt.tight_layout()
plt.savefig('aaf_plot', dpi=200)