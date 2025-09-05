#!/usr/bin/env python3

import sys
import re
import os
import pandas as pd
import subprocess as sp
import glob

# Read input analysis csv and bam files
nextclade_tsv = sys.argv[1]
#bam_file = sys.argv[2]
#sample_name = sys.argv[3]

# Set threshold 
#mix_thresh = sys.argv[4] # This is the proprtion of the major allele
#depth_thresh = 20 # Minimum read depth
#pos_thresh = 4 # Minimum number of positions with mixtures to flag by bammix

meta_analysis_csv = pd.read_csv(nextclade_tsv, sep='\t', usecols=lambda column: column != 'index')

# Get mutations & N's from the meta analysis csv
nextclade_snps = meta_analysis_csv["substitutions"] # mutations
nextclade_Ns = meta_analysis_csv["missing"] # N's

# Drop NaN rows
nextclade_snps.dropna(inplace=True)
nextclade_Ns.dropna(inplace=True)
#print(nextclade_snps)
#print(nextclade_Ns)

# Split nextclade_snps by ',' and remove letters from each element
nextclade_snps = nextclade_snps.str.split(',').apply(
    lambda x: [re.sub('[A-Z]','',i) for i in x]
)

# Split nextclade_Ns by ',' and isolate N's from single positions
nextclade_Ns = nextclade_Ns.str.split(',').apply(
    lambda x: [i for i in x if "-" not in i]
)

# Append all mutations and single N's to a list,
# exclude nan's and empty strings, and sort
all_snps = []
all_single_Ns = []
all_snps_and_single_Ns = []
nextclade_snps = nextclade_snps.apply(lambda x: all_snps.extend(x))
nextclade_Ns = nextclade_Ns.apply(lambda x: all_single_Ns.extend(x))
all_snps = list(set(all_snps))
all_single_Ns = list(set(all_single_Ns))
#print(all_snps)
#print(all_single_Ns)
all_snps_and_single_Ns.extend(all_snps)
all_snps_and_single_Ns.extend(all_single_Ns)
all_snps_and_single_Ns = [
    i for i in all_snps_and_single_Ns if i != '' and i != 'nan'
]
all_snps_and_single_Ns = [int(i) for i in all_snps_and_single_Ns]
all_snps_and_single_Ns = list(set(all_snps_and_single_Ns))
all_snps_and_single_Ns.sort()
#all_snps_and_single_Ns = ' '.join(all_snps_and_single_Ns)
all_snps_and_single_Ns = ' '.join(map(str, all_snps_and_single_Ns))
#all_snps_and_single_Ns = [int(s) for s in all_snps_and_single_Ns_string.split()]
print(all_snps_and_single_Ns)

# convert only for subprocess (command-line args must be str)
#positions = list(map(str, all_snps_and_single_Ns))

#cmd = ["bammix", "-b", bam_file, "-p"] + positions + ["-o", sample_name]
# Run bammix on bam file
#cmd = ["bammix", "-b", bam_file, "-p", all_snps_and_single_Ns, "-o", sample_name]
#print("Running:", " ".join(cmd))
#sp.run(cmd, check=True)
