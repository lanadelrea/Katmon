#!/usr/bin/env python3

import sys
import csv
import subprocess as sp

csv_file = sys.argv[1]
in_dir = sys.argv[2]

with open(csv_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    next(csv_reader)
    sample_names = [row[0].lstrip("./") for row in csv_reader]
    for sample_name in sample_names:
        print(sample_name)
        filtered_bam = f'samtools view -b -q 30 {in_dir}/{sample_name}.bam > {sample_name}_filtered.bam'
        filtered_bam_output = sp.check_output(filtered_bam, shell=True).decode()