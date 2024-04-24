#!/usr/bin/env nextflow

process virstrain {
        tag "Identifying lineage assignment for ${fastqPath.baseName} using Virstrain"
        container 'lanadelrea/virstrain:v0.2.0'

        publishDir (
        path: "${params.out_dir}/03-Virstrain",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        path fastqPath

        script:
        """
        python /usr/src/app/VirStrain.py \
        -i ${params.in_dir}/${fastqPath} \
        -d ${baseDir}/bin/VirStrain/Custom_DB \
        -o ${params.out_dir}/03-Virstrain/${fastqPath.baseName}
        """
}