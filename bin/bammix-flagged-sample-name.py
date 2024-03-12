import sys
import csv

csv_file = sys.argv[1]

with open(csv_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    next(csv_reader)
    sample_names = [row[0].lstrip("./") for row in csv_reader]

with open('samples-flagged.csv', 'w') as txt_file:
    txt_file.write("sample_name\n")
    for sample_name in sample_names:
        txt_file.write(sample_name + '\n')
