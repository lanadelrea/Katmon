#!/bin/bash

# Specify the directory containing the files
directory="/mnt/c/COINFECTION_PAPER/Samples/Batch28_Copy"

# Navigate to the directory
cd "$directory" || exit

# Loop through each file in the directory
for file in *barcode*.primertrimmed.rg.sorted.bam; do
    # Extract the barcode number from the file name
    barcode=$(echo "$file" | grep -oE 'barcode[0-9]+')

    # Rename the .bam file to include only the barcode number
    mv "$file" "${barcode}.bam"

    # Rename the corresponding .bai file
    mv "${file}.bai" "${barcode}.bam.bai"

    # Find and rename the corresponding .fasta file
    fasta_file=$(find . -type f -name "*${barcode}*.consensus.fasta")
    if [ -n "$fasta_file" ]; then
        mv "$fasta_file" "${barcode}.fasta"
    fi

    # Find and rename the corresponding .fastq file
    fastq_file=$(find . -type f -name "*${barcode}*.fastq")
    if [ -n "$fastq_file" ]; then
        mv "$fastq_file" "${barcode}.fastq"
    fi
done
