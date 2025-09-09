// enable dsl2
nextflow.enable.dsl=2

//import modules
include { pangolin } from '../modules/01-lineageAssignment.nf'
include { nextclade } from '../modules/01-lineageAssignment.nf'
include { summary } from '../modules/01-lineageAssignment.nf'

workflow lineage_assignment {
    take: 
        ch_cat_fasta

    main:
        pangolin( ch_cat_fasta )
        nextclade( ch_cat_fasta, params.SC2_dataset )
        summary( pangolin.out.pangolin_csv, nextclade.out.nextclade_tsv )

    emit:
        nextclade_tsv = nextclade.out.nextclade_tsv
        lineage_summary = summary.out.tsv
}

