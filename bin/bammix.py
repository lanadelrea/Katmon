#!/usr/bin/env python3

import sys
import re
import os

import pandas as pd
import subprocess as sp
import glob

# Read input analysis csv, set path to artic analysis folder.
#meta_folder = sys.argv[1]
nextclade_tsv = sys.argv[1]
analysis_folder = sys.argv[2]

# Set threshold
#mix_thresh = 0.8
mix_thresh = sys.argv[3] # This is the proportion of the major allele
depth_thresh = 20 # read depth
pos_thresh = 4 # minimum number of positions with mixtures to flag by bammix

meta_analysis_csv = pd.read_csv(nextclade_tsv, sep='\t', usecols=lambda column: column != 'index')

# Get mutations & N's from the meta analysis csv
#nextclade_snps = meta_analysis_csv.str[15]  # mutations
#nextclade_Ns = meta_analysis_csv.str[29]  # N's
nextclade_snps = meta_analysis_csv["substitutions"]  # mutations
nextclade_Ns = meta_analysis_csv["missing"]  # N's

# Drop NaN rows.
nextclade_snps.dropna(inplace=True)
nextclade_Ns.dropna(inplace=True)
print(nextclade_snps)
print(nextclade_Ns)


# Split nextclade_snps by ',' and remove letters from each element.
nextclade_snps = nextclade_snps.str.split(',').apply(
    lambda x: [re.sub('[A-Z]','',i) for i in x]
    )

# Split nextclade_Ns by ',' and isolate N's from single positions.
nextclade_Ns = nextclade_Ns.str.split(',').apply(
    lambda x: [i for i in x if "-" not in i]
    )


# Append all mutations and single N's to a list,
# exclude nan's and empty strings, and sort.
all_snps = []
all_single_Ns = []
all_snps_and_single_Ns = []
nextclade_snps = nextclade_snps.apply(lambda x: all_snps.extend(x))
nextclade_Ns = nextclade_Ns.apply(lambda x: all_single_Ns.extend(x))
all_snps = list(set(all_snps))
all_single_Ns = list(set(all_single_Ns))
print(all_snps)
print(all_single_Ns)
all_snps_and_single_Ns.extend(all_snps)
all_snps_and_single_Ns.extend(all_single_Ns)
all_snps_and_single_Ns = [
    i for i in all_snps_and_single_Ns if i != '' and i != 'nan'
    ]
all_snps_and_single_Ns = [int(i) for i in all_snps_and_single_Ns]
all_snps_and_single_Ns = list(set(all_snps_and_single_Ns))
all_snps_and_single_Ns.sort()
all_snps_and_single_Ns = [str(i) for i in all_snps_and_single_Ns]
all_snps_and_single_Ns = ' '.join(all_snps_and_single_Ns)


# Glob for names of all bam files in analysis folder.
bam_files = sp.check_output(
    f'find {analysis_folder} -name "*.bam" -type f', shell=True
    ).decode('utf-8').split('\n')
bam_files = [b for b in bam_files if b != '']


# Check if BAM files found
if bam_files:
    print("BAM files found.")
    for bam_file in bam_files:
        print(bam_file)
else:
    print("No BAM files found.")
#except sp.CalledProcessError as e:
#    print(f"Error:{e}")

print(f"Path for finding BAM files: {analysis_folder}") # Look at path to BAM files

##############################################################################################################################
##############################################################################################################################    
##############################################################################################################################

# Glob for names of all BAM files in the bam folder
#bam_files = glob.glob(f"{analysis_folder}/**/*.bam", recursive=True) 

# Call bammix command from shell with all_snps_and_single_Ns as input. We're looking if there are mixture in these positions.
# Use run name, barcode, and central_id as prefix.
for bam in bam_files:
    #prefix  = re.sub('.bam','',bam.split('/')[1])
    #prefix_match = re.search(r'/([^/]+)\.(?:primertrimmed\.rg\.sorted\.bam|.bam)', bam)
    filename = os.path.basename(bam)
    prefix_match= re.match(r'(.+?)\.bam', filename)
    if prefix_match:
        prefix = prefix_match.group(1)
    else:
        print(f"Error: Could not extract prefix from {bam}")
        continue
    print("Prefix:", prefix)
    bammix_cmd = f'bammix -b {bam} -p {all_snps_and_single_Ns} -o {prefix}'
    bammix_output = sp.check_output(bammix_cmd, shell=True).decode()
    print(bammix_output)
    # TO DO: Split list of positions if >N and cannot fit in one plot.
##############################################################################################################################
##############################################################################################################################    
##############################################################################################################################
# Glob for names of all bammix csv files.
bammix_csv_files = sp.check_output(
    f'find ./ -name "*_position_base_counts.csv"', shell=True
    ).decode('utf-8').split('\n')
bammix_csv_files = [i for i in bammix_csv_files if i != '']

if bammix_csv_files:
    print("CSV files found:")
    for csv_file in bammix_csv_files:
        print(csv_file)
else:
    print("No CSV files found.")

# Glob for names of all bammix CSV files
#bammix_csv_files = glob.glob(f"{bammix_process}/**/*_position_base_counts.csv", recursive=True)

# Read all bammix csv files as df's.
bammix_csv = [pd.read_csv(f) for f in bammix_csv_files]


# Iterate over bammix df's, keep track which barcodes had bad mixtures,
# (i.e. >4 positions with >20% mixture at >20x depth),
# at what positions, and what bases.
bad_mixtures = []
bad_mixtures_positions = []
bad_mixtures_base_prop = []
for nam, df in zip(bammix_csv_files, bammix_csv):
    # Get barcode from file name.
    name = re.sub('_position_base_counts.csv','',nam)
    # Get df of nucleotide positions.
    df_pos = df[["Position"]]
    # Get positions with >20% mixture from df.
    df_prop = df[["A_proportion", "C_proportion", "G_proportion", "T_proportion"]]
    # Get largest number for each row in df.
    df_prop_max = df_prop.max(axis=1)
    # Select rows with less than 0.80 value in df_prop_max
    df_prop_max_bool = df_prop_max < mix_thresh
    df_prop_below80max = df_prop[df_prop_max_bool]
    # Flag only if depth is >20 reads.
    df_depth = df[["Total_reads"]]
    df_depth_20x = df_depth["Total_reads"] > depth_thresh
    df_prop_below80max_20xdepth = df_prop_below80max[df_depth_20x]
    second_max_value = df_prop_below80max_20xdepth.apply(lambda row: row.nlargest(2).iloc[-1], axis=1)
    # Flag only if depth is >20 reads, rows have <0.80 value in df_prop_max, and second max value > 0.20.
    df_flagged = df_prop_below80max_20xdepth[(df_prop_max_bool) & df_depth_20x & (second_max_value > 0.20)]
    if df_flagged.empty:
        pass
    elif df_flagged.shape[0] > pos_thresh:
        bad_mixtures.append(name)
        bad_mixtures_positions.append(df_pos.loc[df_flagged.index, "Position"].tolist())
        bad_mixtures_base_prop.append(df_flagged.values.tolist())

# Write table of flagged barcodes, positions, base proportions.
bammix_bad_mixtures_csv = pd.DataFrame()
bammix_bad_mixtures_csv['barcode'] = bad_mixtures
bammix_bad_mixtures_csv['positions'] = bad_mixtures_positions
bammix_bad_mixtures_csv['base_prop'] = bad_mixtures_base_prop
bammix_bad_mixtures_csv.to_csv("bammixFlags.csv", index=False)

# Debugging messages
print("Number of CSV files found:", len(bammix_csv_files))
print("Number of DataFrames in bammix_csv:", len(bammix_csv))

# Concatenate all raw bammix csv files into one csv for reference.
bammix_csv_pos_col = bammix_csv[0]["Position"]
for i in range(len(bammix_csv)):
    bammix_csv[i].columns = [
        f'{bammix_csv_files[i].split("/")[1]}_{c}'
        for c in bammix_csv[i].columns
        ]

bammix_csv = pd.concat(bammix_csv, ignore_index=False, axis=1)
position_cols = [c for c in bammix_csv.columns if 'Position' in c]
bammix_csv.drop(position_cols, axis=1, inplace=True)
bammix_csv = pd.concat([bammix_csv_pos_col,bammix_csv], ignore_index=False, axis=1)
bammix_csv.sort_values("Position", ascending=True, inplace=True, na_position='last')
bammix_csv.to_csv("concat_raw_bammix.csv", index=False)