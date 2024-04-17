#!/usr/bin/env nextflow

process concat {
        tag "Concatenating fasta files into all_sequences.fasta"

        publishDir (
        path: "${params.out_dir}/01-LineageAssignment",
        mode: 'copy',
        overwrite: 'true'
        )

        output:
        path ('*.fasta'), emit: fasta

        script:
        """
        cat ${PWD}/${params.in_dir}/*.fasta > all_sequences.fasta
        """
}

process pangolin {
        cpus 1
        container 'staphb/pangolin:latest'
        tag "Lineage assignment using pangolin tool"

        publishDir (
        path: "${params.out_dir}/01-LineageAssignment",
        mode: 'copy',
        overwrite: 'true',
        )

        input:
        path fasta

        output:
        path ('*.csv'), emit: pangolin_csv

        script:
        """
        pangolin all_sequences.fasta
        """
}

process nextclade {
        cpus 1
        container 'nextstrain/nextclade:latest'
        tag "Lineage assignment using nextclade"

        publishDir (
        path: "${params.out_dir}/01-LineageAssignment",
        mode: 'copy',
        overwrite: 'true',
        )

        input:
        path fasta
        path SC2_dataset

        output:
        path ('*.tsv'), emit: nextclade_tsv

        script:
        """
        nextclade run --input-dataset ${SC2_dataset} --output-tsv=nextclade.tsv ${fasta}
        """
}

process lineage_assignment {
        tag "Create summary table for lineage assignment"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.out_dir}/01-LineageAssignment",
        mode: 'copy',
        overwrite: 'true',
        )

        input:
        path pangolin_csv
        path nextclade_tsv

        output:
        path ('*.tsv'), emit: lineageAssign_tsv

        script:
        """
        python ${baseDir}/bin/lineageAssign_table.py ${pangolin_csv} ${nextclade_tsv}
        """
}