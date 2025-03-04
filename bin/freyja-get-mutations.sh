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

# Concatenate all generated TSVs if they exist
#cat ${sample_name}_*.tsv > ${sample_name}_mutations.tsv
#if ls "${sample_name}"_*.tsv 1> /dev/null 2>&1; then
#    cat "${sample_name}"_*.tsv > "${sample_name}_mutations.tsv"
#    echo "Merged mutations saved to ${sample_name}_mutations.tsv"
#else
#    echo "No TSV files were generated."
#fi

#   command:       freyja-get-lineage.sh ${lineage_list} ${sample} ${annot} ${ref}