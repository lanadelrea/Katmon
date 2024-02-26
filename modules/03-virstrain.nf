#!/usr/bin

process virstrain {
        tag "Identifying viral strains using Virstrain"

        publishDir (
        path: "${params.out_dir}/03-Virstrain",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        path fastqPath

        output:
        tuple val(fastqPath.name), path ('*.tsv'), path ('*.png'), path ('*.html'), emit: virstrain_result
        
        script:
        """
        mkdir -p ${params.out_dir}/03-Virstrain/${fastqPath.baseName}
        
        #!/usr/bin/env python
        python ${baseDir}/bin/VirStrain/VirStrain.py -i ${params.in_dir}/${fastqPath} -d ${baseDir}/bin/VirStrain/Custom_DB -o ${params.out_dir}/03-Virstrain/${fastqPath.baseName}
        """
}