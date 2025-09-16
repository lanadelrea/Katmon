#!/usr/bin/env python3

import sys
import pandas as pd

# Input the collection of sample name and csv path
# Alternating sample name and csv
args = sys.argv[1:]

sample_names = [] # Store names of samples which have positions flagged

# Sanity check: make sure that the arguments come in pairs
if len(args) % 2 != 0:
    raise ValueError("Arguments must come in pairs: <sample> <csv file>")

# Loop through sample/csv pairs
for i in range(0, len(args), 2):
    sample_name = args[i]
    csv_file = args[i + 1]

    try:
        df = pd.read_csv(csv_file)

        if 'positions' not in df.columns: # Check if 'positions' column exists
            print(f"Skipping {csv_file} (no 'positions' column)")
            continue

        # Keep only non-null values in df['positions']
        df_positions = df["positions"].dropna()

        # Check if there are more than 4 rows/postions flagged with mixture
        if len(df_positions) > 4:
            sample_names.append(sample_name)

    except Exception as e:
        print(f"Error processing {csv_file}: {e}")

# Write all flagged sample names into a single file
with open("list_samples_bammix_flagged.txt", "w") as out:
    for name in sample_names:
        out.write(name + "\n")