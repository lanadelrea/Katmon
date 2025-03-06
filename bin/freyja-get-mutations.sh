#!/bin/bash

list_of_lineage="$1"
sample_name="$2"
annot="$3"
ref="$4"

i=2

while [ ${i} -lt 4 ]; do
        # Get the lineage from line 2 of lineage list
        lineage=$(sed -n "${i}p" "${list_of_lineage}")
        echo "${lineage}"

        # Ensure that lineage is not empty
        if [ -n "${lineage}" ]; then
            freyja get-lineage-def ${lineage} --annot ${annot} --ref ${ref} --output ${sample_name}_${lineage}.tsv
        else
            echo "No lineage found at line ${i}, exiting."
            break
        fi

        # Add one to get line 3 from the list
        i=$((i + 1))
done

#   command:       freyja-get-lineage.sh ${lineage_list} ${sample} ${annot} ${ref}
# Writtten by Adeliza Realingo (04 March 2025)