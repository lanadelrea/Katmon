#!/usr/bin/env nextflow

process virstrain {
        tag "Identifying lineage assignment for ${fastqPath.BaseName} using VirStrain"
        container 'lanadelrea/virstrain:v0.3.0'

        publishDir (
        path: "${params.out_dir}/03-VirStrain/VirStrain_txt",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        path(virstrain_database)
        path(fastqPath)

        output:
        path('*.txt'), emit: txt

        script:
        """
        virstrain \
        -i ${fastqPath} \
        -d ${virstrain_database}

        mv ./VirStrain_Out/VirStrain_report.txt ${fastqPath.BaseName}.txt
        """
}

process virstrain_summary {
        tag "Summarizing VirStrain results"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.out_dir}/03-VirStrain",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (virstrain_txt_dir)
        path (txt_files)

        output:
        path ('virstrainSummary.tsv'), emit: tsv

        script:
        """
        virstrain_table.py ${virstrain_txt_dir}
        """
}