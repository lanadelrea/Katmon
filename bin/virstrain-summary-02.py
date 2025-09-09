#!/usr/bin/env python3

# by ALR/lanadelrea/09September2025
# Make a summary table of the VirStrain results

import sys
import csv
import os
import pandas as pd


# Input the list of the tsv files collected from the previous process
tsv_files = sys.argv[1:]

results = [] # To store the results

# Append the results to create one table
for tsv_file in tsv_files:
    try:
        df = pd.read_csv(tsv_file, sep = '\t')
        if not df.empty:
            virstrain_result = df.iloc[0] # Takes the first row containing the virstrain result
            results.append(virstrain_result)
    except Exception as e:
        print(f"Error processing {tsv_file}: {e}")

if results:
    result_df = pd.DataFrame(results)
    result_df.to_csv("VirStrain_Results.tsv", sep='\t', index=False)
    print("Results are summarized in VirStrain_Results.tsv")
else:
    print("No data collected from input files")