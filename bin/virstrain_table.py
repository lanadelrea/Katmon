#!/usr/bin/env python3

import os
import re
import pandas as pd
import sys

# Directory containing the text files
directory_path = sys.argv[1]

# Regular expression pattern to match the first two strings separated by "_"
pattern = r'>([^_]+)_([^_]+)'

# Regular expression pattern to match the strain name in the ">>Other possible strains:" section
strain_pattern = r'>([^_]+)'

# Function to extract the specific substrings from the given line
def extract_substrings(line):
    match = re.search(pattern, line)
    if match:
        return match.group(1), match.group(2)
    return None, None

# Function to extract the strain name from the given line
def extract_strain(line):
    match = re.search(strain_pattern, line.strip())
    if match:
        return match.group(1)
    return None

# Initialize a list to hold data for the DataFrame
data = []

# Function to process each file
def process_file(file_path):
    most_possible_strains = []
    other_possible_strains = []
    capture_most_possible = False
    capture_other_possible = False

    with open(file_path, 'r') as file:
        for line in file:
            stripped_line = line.strip()

            # Check if we are at the start of the ">>Most possible strains:" section
            if stripped_line == ">>Most possible strains:":
                capture_most_possible = True
                capture_other_possible = False
                continue

            # Check if we are at the start of the ">>Other possible strains:" section
            if stripped_line == ">>Other possible strains:":
                capture_other_possible = True
                capture_most_possible = False
                first_capture_done = False
                continue

            # If we are capturing lines and encounter a new section, stop capturing
            if stripped_line.startswith(">>") and not stripped_line.startswith(">>>"):
                capture_most_possible = False
            #if stripped_line.startswith(">"):
            #    capture_other_possible = False

            # Capture the lines in the ">>Most possible strains:" section that start with '>'
            if capture_most_possible and stripped_line.startswith('>'):
                most_possible_strains.append(stripped_line)

            if capture_other_possible:
                if stripped_line == "Can not detect other strains.":
                    other_possible_strains.append(stripped_line)
                    first_capture_done = True
                elif stripped_line.startswith('>') and not stripped_line.startswith('>>') and not first_capture_done:
                    other_possible_strains.append(stripped_line)
                    first_capture_done = True

    # Extract the desired substrings from each captured line in ">>Most possible strains:" section
    most_possible_strains_extracted = []
    for strain in most_possible_strains:
        first_str, second_str = extract_substrings(strain)
        if first_str and second_str:
            most_possible_strains_extracted.append(f"{first_str}_{second_str}")

    # Extract the desired information from each captured line in ">>Other possible strains:" section
    other_possible_strains_extracted = []
    for strain in other_possible_strains:
        if strain == "Can not detect other strains.":
            other_possible_strains_extracted.append(strain)
        else:
            strain_name = extract_strain(strain)
            if strain_name:
                other_possible_strains_extracted.append(strain_name)

    # Combine extracted strains into a single string for each section
    most_possible_strains_combined = "; ".join(most_possible_strains_extracted)
    other_possible_strains_combined = "; ".join(other_possible_strains_extracted)

    # Add the data to the list
    file_name = os.path.basename(file_path)
    sample_name = file_name.split('_')[0]
    data.append([sample_name, most_possible_strains_combined, other_possible_strains_combined])

# Process each file in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".txt"):  # Assuming the files are .txt files
        file_path = os.path.join(directory_path, filename)
        process_file(file_path)
 
# Create a DataFrame from the data
df = pd.DataFrame(data, columns=["Sample", "Most Possible Strain", "Other Possible Strain"])

df.to_csv("virstrainSummary.tsv", sep='\t', index=False)