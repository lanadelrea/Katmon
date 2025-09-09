#!/usr/bin/env nextflow

process virstrain_process {
        tag "Identifying lineage assignment for ${sample} using VirStrain"
        container 'lanadelrea/virstrain:v0.3.0'

        publishDir (
        path: "${params.out_dir}/03-VirStrain/${sample}",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        path(virstrain_database)
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path ("VirStrain_Out/*.txt"), emit: txt
        path ("VirStrain_Out/*.html")
        path ("VirStrain_Out/*.csv")

        script:
        """
        virstrain \
        -i ${fastq} \
        -d ${virstrain_database}
        """
}

process wrangling {
        tag "Data wrangling of information from VirStrain output text"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.out_dir}/03-VirStrain/${sample}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(virstrain_txt)

        output:
        path ('*.tsv'), emit: tsv

        script:
        """
        virstrain-wrangling-01.py ${sample} ${virstrain_txt}
        """
}

process summary {
        tag "Creating a summary table of VirStrain results"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.out_dir}/03-VirStrain",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (tsv)
        
        output:
        path ('*.tsv'), emit: tsv

        script:
        """
        virstrain-summary-02.py ${tsv}
        """
}