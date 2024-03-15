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

# Output files
read_depth_path = sys.argv[5]   # 'read_depths/PH-RITM-1395.tsv'
plot_path = sys.argv[6]          # 'PH-RITM-1395-final.png'
###########################################################################

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