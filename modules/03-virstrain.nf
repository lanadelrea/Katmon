#!/usr/bin/env nextflow

process virstrain {
        tag "Identifying lineage assignment for ${fastqPath.baseName} using VirsStrain"
        container 'lanadelrea/virstrain:latest'

        publishDir (
        path: "${params.out_dir}/03-VirStrain",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        path fastqPath

        output:
        tuple val(fastqPath.baseName), path ('*.txt'), emit: virstrain_txt

        script:
        """
        virstrain \
        -i ${fastqPath} \
        -d $PWD/CoPi/assets/Custom_DB \
        -o $PWD/${params.out_dir}/03-Virstrain/${fastqPath.baseName}
        """
}

process virstrain_summary_txt {
        tag "Summarizing VirStrain results"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.out_dir}/03-VirStrain",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(txt)

        output:
        path ('*.tsv'), emit: virstrain_tsv

        script:
        """   
        mkdir -p $PWD/${params.out_dir}/03-Virstrain/VirStrain_txt
        cp  $PWD/${params.out_dir}/03-Virstrain/${sample}/*.txt ${sample}_VirStrain.txt
        mv  ${sample}_VirStrain.txt $PWD/${params.out_dir}/03-Virstrain/VirStrain_txt

        virstrain_table.py $PWD/${params.out_dir}/03-Virstrain/VirStrain_txt
        """
}