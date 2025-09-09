// enable dsl2
nextflow.enable.dsl=2

// import modules
include { virstrain_process } from '../modules/03-virstrain.nf'
include { wrangling } from '../modules/03-virstrain.nf'
include { summary } from '../modules/03-virstrain.nf'

workflow virstrain {
    take:
        ch_fastq_files
    
    main:
        virstrain_process(params.virstrain_database, ch_fastq_files)
        wrangling(virstrain_process.out.txt)
        summary(wrangling.out.tsv.collect())

    emit:
        virstrain_summary = summary.out.tsv
}