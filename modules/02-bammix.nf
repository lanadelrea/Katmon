#!/usr/bin/env nextflow

process bammix {
        tag "Looking for positions with nucleotide mixtures using bammix"
        container 'ufuomababatunde/bammix:v1.0.0'

        publishDir "${params.out_dir}/${task.process.replaceAll(":","_")}", pattern: "*.csv", mode: 'copy'
        publishDir "${params.out_dir}/${task.process.replaceAll(":","_")}", pattern: "*.pdf", mode: 'copy'

        input:
        path(nextclade_tsv)
        path(bam)
        path(bam_index)

        output:
        path '*.csv'
        path '*.pdf'

        script:
        """
        bammix.py ${nextclade_tsv} ${params.in_dir}
        """
}