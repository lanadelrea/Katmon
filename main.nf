// enable dsl2
nextflow.enable.dsl=2

//import modules
include { concat } from './modules/01-lineageAssignment.nf'
include { pangolin } from './modules/01-lineageAssignment.nf'
include { nextclade } from './modules/01-lineageAssignment.nf'
include { bammix } from './modules/02-bammix.nf'
include { bam_filter } from './modules/02-bammix.nf'
include { virstrain } from './modules/03-virstrain.nf'
include { freyja } from './modules/04-freyja.nf'
include { freyja_demix } from './modules/04-freyja.nf'
include { freyja_aggregate } from './modules/04-freyja.nf'
include { freyja_plot } from './modules/04-freyja.nf'
include { makevcf } from './modules/05-makeVCF.nf'
include { bammixplot } from './modules/06-plots.nf'
include { aafplot_mutations } from './modules/06-plots.nf'
include { aafplot_amplicons } from './modules/06-plots.nf'
//include { ampliconsorting } from './modules/06-ampliconSorting.nf'
//include { report } from '.modules/07-report.nf'

workflow {
        ch_bam_file = Channel
                    .fromPath("${params.in_dir}/**.bam", type: 'file')
                    .ifEmpty { error "Cannot find any BAM files on ${params.in_dir}"}
        ch_bam_index = Channel
                    .fromPath("${params.in_dir}/**.bai", type: 'file')
                    .ifEmpty { error "Cannot find any BAM index files on ${params.in_dir}"}
        ch_fastq = Channel
                    .fromPath("${params.in_dir}/**.fastq", type: 'file')
                    .ifEmpty { error "Cannot find any fastq files" on ${params.in_dir}}
        ch_fasta = Channel
                    .fromPath("${params.in_dir}/**.fasta", type: 'file')
                    .ifEmpty { error "Cannot find any fasta files" on ${params.in_dir}}

        main:
               ch_bam_file.map { bamfilePath -> tuple(bamfilePath) }
               ch_fastq.map { fastqPath -> tuple(fastqPath) }
               ch_fasta.map { fastaPath -> tuple(fastaPath) }

               concat()
               pangolin( concat.out.fasta )
               nextclade( concat.out.fasta )
               bammix ( nextclade.out.nextclade_tsv, ch_bam_file, ch_bam_index )
               bam_filter ( bammix.out.bammixflagged_csv)
               virstrain ( ch_fastq )
               freyja( ch_bam_file )
               freyja_demix( freyja.out.freyja_variants )
               freyja_aggregate( freyja_demix.out.tsv_demix.collect().view() )
               freyja_plot( freyja_aggregate.out.freyja_aggregated_file )
               makevcf( bam_filter.out.filtered_bam )
               bammixplot (makevcf.out.filtered_vcf.collect())
               aafplot_mutations( makevcf.out.filtered_vcf.collect())
               aafplot_amplicons( aafplot_mutations.out.aafplot_mut.collect()) 
}