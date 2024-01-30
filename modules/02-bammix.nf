#!/usr/bin/env nextflow

process bammix {
        tag "Looking for positions with nucleotide mixtures using bammix"
        container 'ufuomababatunde/bammix:v1.0.0'

        publishDir (
        path: "${params.out_dir}/02-BammixFlags",
        pattern: "*.csv",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path(nextclade_tsv)
        path(bam)
        path(bam_index)

        output:
        path '*.csv'
        path "flagged_barcode_positions_proportions.csv", emit: bammix_summary

        script:
        """
        bammix.py ${nextclade_tsv} ${params.in_dir}
        """
}
