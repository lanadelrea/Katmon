#!/usr/bin/env python3

import csv
import argparse
import ast

def extract_bammix(int_list, output_file):
    # Convert eaxh integer into a single-element iterable
    transposed_list = [[x] for x in int_list]

    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(transposed_list)

# Create an argument parser
parser = argparse.ArgumentParser(description='Extract a list of integers and output to a transposed CSV file')

# Add the list argument
parser.add_argument('--list', type=str, help='Input list of values')

# Add the output file argument
parser.add_argument('--out', type=str, help='Output CSV file')

# Parse the command-line arguments
args = parser.parse_args()

# Extract the list from the string argument
input_list = ast.literal_eval(args.list)

# Example usage
output_csv_file = args.out
extract_bammix(input_list, output_csv_file)
print("Extraction complete. Please check the transposed output file:", output_csv_file)