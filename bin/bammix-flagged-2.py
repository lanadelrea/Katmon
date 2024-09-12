#!/usr/bin/env python3

import sys
import csv
import subprocess as sp
import os

# Input the CSV file of the bammix flagged samples and the input directory
csv_file = sys.argv[1]

# Reads the CSV file
with open(csv_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    next(csv_reader)
    sample_names = [row[0].lstrip("./") for row in csv_reader]  # Gets the sample names flagged by bammix
    
    for sample_name in sample_names:
        print(f'{sample_name}')