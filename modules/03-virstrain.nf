#!/usr/bin

process virstrain {
        tag "Identifying viral strains using Virstrain"

        publishDir (
        path: "${params.out_dir}/03-Virstrain",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        path(fasta)

        output:
        path '*.tsv'
        path '*.png'
        path '*.html'
        
        script:
        """
        python Virst
        python Virstrain_contig.py -i ${fasta} -v 1 -d ${Virstrain_DB} -o ${params.out_dir}
        """
}