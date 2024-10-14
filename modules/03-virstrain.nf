#!/usr/bin/env nextflow

process virstrain {
        tag "Identifying lineage assignment for ${fastqPath.SimpleName} using VirsStrain"
        container 'lanadelrea/virstrain:v0.3.0'

        publishDir (
        path: "${params.out_dir}/03-VirStrain/",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        path(fastqPath)
//        path(fasta)

        output:
        path('*.txt'), emit: txt

        script:
        """
        virstrain \
        -i ${fastqPath} \
        -d $PWD/Katmon/assets/Custom_DB \
        -o ${params.out_dir}/03-Virstrain/${fastqPath.SimpleName}

        cp ${params.out_dir}/03-Virstrain/${fastqPath.SimpleName}/VirStrain_report.txt ${fastqPath.SimpleName}_VirStrain.txt
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
        path (txt_files)

        output:
        path ('*.tsv'), emit: tsv

        script:
        """
        mkdir VirStrain_txt
        cp -r ${txt_files} VirStrain_txt/
        virstrain_table.py VirStrain_txt/
        """
}