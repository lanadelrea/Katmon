#!/usr/bin/env python

import argparse

from ast import literal_eval
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pysam
import sys

cmap=sns.color_palette('pastel')

######################## F I R S T  P L O T ##############################
# Input files
vcf_path = sys.argv[1]          # 'vcf/PH-RITM-1395.vcf.gz'
mutation_path = sys.argv[2]     # 'mutations-final-update.tsv'
bam_file = sys.argv[3]          # 'PH-RITM-1395.primertrimmed.rg.sorted.bam'

sample = sys.argv[4]            # 'PH-RITM-1395'

# Output
read_depth_path = sys.argv[5]   # 'read_depths/PH-RITM-1395.tsv'
plot_path = sys.argv[6]          # 'PH-RITM-1395-final.png'
#aaf_plot = sys.argv[7]          # 'PH-RITM-1395-AAF.png'


def get_lineage(delta_annotation: str, omicron_annotation: str) -> str:
    if delta_annotation == omicron_annotation == None:
        return 'None'
    if delta_annotation == omicron_annotation == 'yes':
        return 'Both'
    else:
        if isinstance(delta_annotation, str) and delta_annotation.lower() == 'yes':
            return 'Delta'
        if isinstance(omicron_annotation, str) and omicron_annotation.lower() == 'yes':
            return 'Omicron'
        
vcf_file = pd.read_csv(vcf_path, header=None, sep='\t', comment='#')

# Extract values in DP4 in column 7 of the vcf file
ref_forward = vcf_file[7].str.extract(r'DP4=(\d+)')
ref_forward = ref_forward.replace(np.nan, 0)
ref_forward = ref_forward.astype(int)

ref_reverse = vcf_file[7].str.extract(r'DP4=\d+(\d+)')
ref_reverse = ref_reverse.replace(np.nan, 0)
ref_reverse = ref_reverse.astype(int)

alt_forward = vcf_file[7].str.extract(r'DP4=\d+,\d+,(\d+),\d+')
alt_forward = alt_forward.replace(np.nan, 0)
alt_forward = alt_forward.astype(int)

alt_reverse = vcf_file[7].str.extract(r'DP4=\d+,\d+,\d+,(\d+)')
alt_reverse = ref_reverse.replace(np.nan, 0)
alt_reverse = alt_reverse.astype(int)

vcf = vcf_file

vcf['TotalReads'] = ref_forward + ref_reverse + alt_forward + alt_reverse
vcf['BaseCalledReadsWithVariant'] = alt_forward + alt_reverse

vcf.rename(columns={ 0: 'chrom', 1: 'pos', 2: 'id', 3: 'ref', 4: 'alt', 5: 'qual', 6: 'filter',
       7 : 'info', 8: 'format', 9 : 'sample', 'TotalReads' : 'dp', 'BaseCalledReadsWithVariant' : 'ad'} , inplace=True)
vcf.set_index(['pos', 'ref', 'alt'], inplace=True)
vcf = vcf[~vcf.index.duplicated(keep='first')]

# Read in curated mutations 'mutations.tsv'
mapping = pd.read_csv(mutation_path, sep='\t')
mapping.dropna(axis='index', how='all', subset=['Delta', 'Omicron'], inplace=True)

# Data wrangling
mapping['VCF'] = mapping['VCF'].apply(literal_eval)
mapping['VCF'] = mapping['VCF'].apply(lambda x: (int(x[0]), x[1], x[2]))
mapping['position'] = mapping['VCF'].apply(lambda x: x[0])
mapping['ref'] = mapping['VCF'].apply(lambda x: x[1])
mapping['alt'] = mapping['VCF'].apply(lambda x: x[2])

# Adding annotations
mapping['aa'] = mapping.apply(lambda x: f'{x.gene}:{x.protein_change}', axis=1)
mapping['lineage'] = mapping.apply(lambda x: get_lineage(x.Delta, x.Omicron), axis=1)
mapping = mapping.set_index(['position', 'ref', 'alt'])

# Loop through curated mutations, saving read depth information to file.
with open(read_depth_path, 'w') as fout:
    # Print column headers.
    print('name', 'variant', 'aa', 'alt_ad', 'dp', 'aaf', 'mutation_label', sep='\t', file=fout)
     # Loop through curated mutations, collecting read depth information.
    for mut in mapping.index.tolist():
        if mut in vcf.index:
            # If position and mutation match, print the allele depths
            # alt_ad = list(map(int, vcf.at[mut,'ad'].tolist()))
            alt_ad = vcf.loc[mut, 'ad']
            dp = vcf.loc[mut, 'dp']
            print(sample,
                  mut,
                  mapping.loc[mut,'aa'].values[0],
                  alt_ad,
                  dp,
                  alt_ad/dp,
                  mapping.loc[mut,'lineage'].values[0],
                  sep='\t',
                  file=fout)
        elif mut not in vcf.index:
            # This is a reference call, so there is not AD field, only the DP field.
            print(sample,
                  mut,
                  mapping.loc[mut,'aa'].values[0],
                  0,
                  0,
                  0,
                  mapping.loc[mut,'lineage'].values[0],
                  sep='\t',
                  file=fout)
        else:
            # Otherwise, find the record matching the position of the curated mutation.
            for record in vcf.index.tolist():
                if record == (mut[0], mut[1], '.'):
                    # This is a reference call, so there is not AD field, only the DP field.
                    print(sample,
                          mut,
                          mapping.loc[mut,'aa'].values[0],
                          0,
                          int(vcf.at[record,'dp']),
                          0,
                          mapping.loc[mut,'lineage'].values[0],
                          sep='\t',
                          file=fout)
                elif (record[0] == mut[0]) and (record[1] == mut[1]):
                    # Not a reference call, but the alternate allele is different, or there are multiple alternate alleles.
                    ads = vcf.at[record,'ad']
                    # Loop through all possible alleles (usually there's 1, but there can sometimes be more).
                    for i, alt in enumerate(record[2].split(',')):
                        new_record = (record[0], record[1], alt)
                        # If one of the alt alleles matches a curated mutation, report the read depths.
                        if new_record == mut:
                            print(sample,
                                  mut,
                                  mapping.loc[mut,'aa'].values[0],
                                  ads[i+1],
                                  sum(ads),
                                  ads[i+1]/sum(ads),
                                  mapping.loc[mut,'lineage'].values[0],
                                  sep='\t',
                                  file=fout)


df = pd.read_csv(read_depth_path, sep='\t')
df_stat = df[df['aaf'] != 0 ]
plt.figure(figsize=(18,4))

color_dict = { 'Both': '#967bb6',
               'Omicron': '#efbab4',
               'Delta': '#1484b4'}
sns.stripplot(x='aa', y='aaf', s=10, hue='mutation_label', palette=color_dict,
              data=df_stat)
plt.ylabel('Allele fraction')
plt.xlabel('')
plt.title(sample)
plt.xticks(rotation=90)
plt.ylim([0,1.1])
plt.legend(bbox_to_anchor=[1,1], ncol=1)
plt.tight_layout()
plt.savefig(plot_path, dpi=200) 

######################## S E C O N D  P L O T ############################

# Input files
primer_scheme_file = sys.argv[7]          # SARS-CoV-2.scheme.bed
primer_scheme = pd.read_csv(primer_scheme_file, sep='\t')

# Output files
aaf_plot = sys.argv[8]          # 'AAFplot_amplicons.png'

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

alt_plot = df

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