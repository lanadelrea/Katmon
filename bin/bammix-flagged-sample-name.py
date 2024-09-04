#!/usr/bin/env python3

import sys
import csv
import subprocess as sp
import os

# Input the CSV file of the bammix flagged samples and the input directory
csv_file = sys.argv[1]
in_dir = sys.argv[2]

# Convert to absolute path if necessary
if not os.path.isabs(in_dir):
    in_dir = os.path.abspath(in_dir)

# Reads the CSV file
with open(csv_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    next(csv_reader)
    sample_names = [row[0].lstrip("./") for row in csv_reader]  # Gets the sample names flagged by bammix
    
    for sample_name in sample_names:
        print(f'Processing {sample_name}')
        
        # Search for the .bam file in the input directory
        bam_file_path = None
        for root, dirs, files in os.walk(in_dir):
            for file in files:
                if file == f'{sample_name}.bam':
                    bam_file_path = os.path.join(root, file)
                    break
            if bam_file_path:
                break

        if bam_file_path:
            filtered_bam_command = f'samtools view -b -q 20 {bam_file_path} > {sample_name}.filtered.bam'  # Filters the bam file of reads with at least Q 12
            sp.run(filtered_bam_command, shell=True)
            print(f'Finished {sample_name}.filtered.bam')
        else:
            print(f'BAM file for {sample_name} not found in {in_dir}')
