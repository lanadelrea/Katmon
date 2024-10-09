// enable dsl2
nextflow.enable.dsl=2

//import modules
include { pangolin } from '../modules/01-lineageAssignment.nf'
include { nextclade } from '../modules/01-lineageAssignment.nf'
include { lineage_assignment } from '../modules/01-lineageAssignment.nf'
include { bammix } from '../modules/02-bammix.nf'
include { bam_filter } from '../modules/02-bammix.nf'
include { makevcf } from '../modules/02-bammix.nf'
include { makevcf_2 } from '../modules/02-bammix.nf'
include { bcftools } from '../modules/02-bammix.nf'
include { virstrain } from '../modules/03-virstrain.nf'
include { virstrain_summary } from '../modules/03-virstrain.nf'
include { freyja } from '../modules/04-freyja.nf'
include { freyja_demix } from '../modules/04-freyja.nf'
include { freyja_aggregate } from '../modules/04-freyja.nf'
include { freyja_plot } from '../modules/04-freyja.nf'
include { bammixplot } from '../modules/05-plots.nf'
include { aafplot_mutations } from '../modules/05-plots.nf'
include { aafplot_mutations_2 } from '../modules/05-plots.nf'
include { aafplot_amplicons } from '../modules/05-plots.nf'
include { ampliconsorting_DeltaReads } from '../modules/06-ampliconSorting.nf'
include { ampliconsorting_OmicronReads } from '../modules/06-ampliconSorting.nf'
include { ampliconsorting_samtools } from '../modules/06-ampliconSorting.nf'
include { ampliconsorting_bgzip } from '../modules/06-ampliconSorting.nf'
include { ampliconsorting_fasta } from '../modules/06-ampliconSorting.nf'
include { ampliconsorting_lineageAssignment_Pangolin } from '../modules/06-ampliconSorting.nf'
include { ampliconsorting_lineageAssignment_Nextclade } from '../modules/06-ampliconSorting.nf'
include { report } from '../modules/07-report.nf'
include { report_no_flag } from '../modules/07-report.nf'
workflow no_bammix_flag {
    ch_bam_file = Channel
                 .fromPath("${params.in_dir}/**.bam", type: 'file')
                 .ifEmpty { error "Cannot find any BAM files on ${params.in_dir}"}
    ch_bam_index = Channel
                 .fromPath("${params.in_dir}/**.bai", type: 'file')
                 .ifEmpty { error "Cannot find any BAM index files on ${params.in_dir}"}
    ch_fastq = Channel
                 .fromPath("${params.in_dir}/**.fastq{.gz,}", type: 'file')
                 .ifEmpty { error "Cannot find any fastq files on ${params.in_dir}"}
    ch_fasta = Channel
                 .fromPath("${params.in_dir}/**.fasta", type: 'file')
                 .ifEmpty { error "Cannot find any fasta files on ${params.in_dir}"}
    main: 
        ch_bam_file.map { bamfilePath -> tuple(bamfilePath) }
        ch_fastq.map { fastqPath -> tuple(fastqPath) }
        ch_cat_fasta = ch_fasta.collectFile(name: 'all_sequences.fasta', newLine: true )

        pangolin( ch_cat_fasta )
        nextclade( ch_cat_fasta, params.SC2_dataset )
        lineage_assignment( pangolin.out.pangolin_csv, nextclade.out.nextclade_tsv )

        freyja(ch_bam_file)

        freyja_demix(freyja.out.freyja_variants)
        freyja_aggregate(freyja_demix.out.tsv_demix.collect())
        freyja_plot(freyja_aggregate.out.freyja_aggregated_file)

        virstrain(ch_fastq, lineage_assignment.out.lineageAssign_tsv)
        ch_virstrain_txt = Channel.fromPath("${params.out_dir}/03-VirStrain", type: 'dir')
        virstrain_summary(ch_virstrain_txt.collect())

        report_no_flag(
            lineage_assignment.out.lineageAssign_tsv,
            freyja_plot.out.freyja_plot,
            virstrain_summary.out.tsv,
            params.report_rmd_no_bammix)
}