// enable dsl2
nextflow.enable.dsl=2

//import modules

include { concat } from './modules/01-lineageAssignment.nf'
include { pangolin } from './modules/01-lineageAssignment.nf'
include { nextclade } from './modules/01-lineageAssignment.nf'
include { bammix } from './modules/02-bammix.nf'
include { virstrain } from './modules/03-virstrain.nf'
//include { freyja } from './modules/04-freyja.nf'
//include { aafplot } from '.modules/05-AAFplot.nf'
//include { ampliconsorting } from './modules/06-ampliconSorting.nf'
//include { report } from '.modules/07-report.nf'

workflow {
        ch_bam_file = Channel
                    .fromPath("${params.in_dir}/**.bam", type: 'file')
                    .ifEmpty { error "Cannot find any BAM files on ${params.in_dir}"}
                    .view()
        ch_bam_index = Channel
                    .fromPath("${params.in_dir}/**.bai", type: 'file')
                    .ifEmpty { error "Cannot find any BAM index files on ${params.in_dir}"}
                    .view()
        ch_fastq = Channel
                    .fromPath("${params.in_dir}/**.fastq", type: 'file')
                    .ifEmpty { error "Cannot find any fastq files" on ${params.in_dir}}
                    .view()

        main:
              ch_bam_file.map { bamfilePath -> tuple(bamfilePath) }
               ch_fastq.map { fastqPath -> tuple(fastqPath) }

               concat()
               pangolin( concat.out.fasta )
               nextclade( concat.out.fasta )
               bammix ( nextclade.out.nextclade_tsv, ch_bam_file, ch_bam_index )
               virstrain ( ch_fastq )
//               freyja( )
}