#!/usr/bin

process virstrain {
        tag "Identifying viral strains using Virstrain"

        publishDir (
        path: "${params.out_dir}/03-Virstrain",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        tuple val(sample), path(fasta)

        output:
        path '*.tsv'
        path '*.png'
        path '*.html'
        
        script:
        """
        python Virstrain_contig.py -i ${fastq} -v 1 -d ./assets/VirStrain/Custom_DB -o ${params.out_dir}/03-Virstrain/${sample}
        """
}