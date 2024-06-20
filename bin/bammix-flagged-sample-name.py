#!/usr/bin/env python3

import sys
import csv
import subprocess as sp

# Input the csv file of the bammix flagged samples and the input directory
csv_file = sys.argv[1]
in_dir = sys.argv[2]

# Reads the csv file 
with open(csv_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    next(csv_reader)
    sample_names = [row[0].lstrip("./") for row in csv_reader] # Gets the sample names flagged by bammix
    for sample_name in sample_names:
        print(sample_name)
        filtered_bam = f'samtools view -b -q 12 {in_dir}/{sample_name}.bam > {sample_name}.bam' # Filters the bam file of reads with at least Q 12 
        filtered_bam_output = sp.check_output(filtered_bam, shell=True).decode()