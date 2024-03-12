#!/usr/bin/env nextflow

process bammix {
        tag "Looking for positions with nucleotide mixtures using bammix"
        container 'ufuomababatunde/bammix:v1.0.0'

        publishDir (
        path: "${params.out_dir}/02-BammixForNucleotideMixtures",
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

process bammix_flagged {
        tag "Getting positions flagged by bammix"
        
        publishDir (
        path: "${params.out_dir}/02-BammixForNucleotideMixtures",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path bammixflagged_csv

        output:
        path 'samples-flagged.csv', emit: bammixflagged_samples

        script:
        """
        python bammix-flagged-sample-name.py ${bammixflagged_csv}
        """
}