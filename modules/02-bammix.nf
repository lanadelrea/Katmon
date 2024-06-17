#!/usr/bin/env nextflow

process bammix {
        tag "Looking for positions with nucleotide mixtures using bammix"
        container 'ufuomababatunde/bammix:v1.1.0'

        cpus 1

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path(nextclade_tsv)
        path(bam)
        path(bam_index)

        output:
        path '*.csv'
        path '*.pdf'
        path 'flagged_barcodes_positions_proportions.csv', emit: bammixflagged_csv

        script:
        """
        bammix.py ${nextclade_tsv} ${params.in_dir}
        """
}

process bam_filter {
        tag "Getting positions flagged by bammix"
        container 'pegi3s/samtools_bcftools:latest'

        publishDir (
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path bammixflagged_csv

        output:
        path ('*.bam'), emit: filtered_bam

        script:
        """
        bammix-flagged-sample-name.py ${bammixflagged_csv} ${params.in_dir}
        """
}