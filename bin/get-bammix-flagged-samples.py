#!/usr/bin/env python3

# by ALR/lanadelrea/03September2025

import sys
import csv
import os
import pandas as pd

# Input the list of the csv files collected from the previous process
csv_files = sys.argv[1:]

sample_names = [] # Store names of samples which have positions flagged

# Check to see if there are positions flagged in each sample
for csv_file in csv_files:
    try:
        df = pd.read_csv(csv_file)

        if 'positions' not in df.columns: # Check to see if the positions column exist (sanity check)
            print(f"Skipping {csv_files} (no 'positions' column)")
            continue

        sample_name = (os.path.basename(csv_file)).split("_")[0]

        # Check in 'positions' column has values, then add to list
        if df['positions'].notnull().any():
            sample_names.append(sample_name)
    except Exception as e:
        print(f"Error processing {csv_file}: {e}")

# Writes the sample names with flagged positions to a file
with open("list_samples_bammix_flagged.txt", "w") as out:
    for name in sample_names:
        out.write(name + "\n")