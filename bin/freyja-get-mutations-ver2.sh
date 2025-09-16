#!/bin/bash

list_of_lineage="$1"
sample_name="$2"
annot="$3"
ref="$4"

count=0   # how many valid lineages we’ve collected so far
i=1       # start from line 1

while [ ${count} -lt 3 ]; do
    lineage=$(sed -n "${i}p" "${list_of_lineage}")

    # If we’ve reached the end of the file, break
    if [ -z "${lineage}" ]; then
        echo "Reached end of list, only found ${count} valid lineage(s)."
        break
    fi

    # Skip lineages with '-like'
    if [[ "${lineage}" == *"-like"* ]]; then
        echo "Skipping ${lineage} (contains -like)"
        i=$((i + 1))
        continue
    fi

    # Process lineage
    echo "Processing ${lineage}"
    freyja get-lineage-def ${lineage} --annot ${annot} --ref ${ref} --output ${sample_name}_${lineage}.tsv

    count=$((count + 1))  # increment valid lineage counter
    i=$((i + 1))          # move to next line
done