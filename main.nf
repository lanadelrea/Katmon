// enable dsl2
nextflow.enable.dsl=2

//import modules
include { pangolin } from './modules/01-lineageAssignment.nf'
include { nextclade } from './modules/01-lineageAssignment.nf'
include { lineage_assignment } from './modules/01-lineageAssignment.nf'

include { bammix } from './modules/02-bammix.nf'
include { bam_filter } from './modules/02-bammix.nf'
include { makevcf } from './modules/02-bammix.nf'
include { makevcf_2 } from './modules/02-bammix.nf'
include { bcftools } from './modules/02-bammix.nf'

include { virstrain } from './modules/03-virstrain.nf'
include { virstrain_summary } from './modules/03-virstrain.nf'

include { freyja } from './modules/04-freyja.nf'
include { freyja_demix } from './modules/04-freyja.nf'
include { freyja_aggregate } from './modules/04-freyja.nf'
include { freyja_plot } from './modules/04-freyja.nf'

include { bammixplot } from './modules/05-plots.nf'
include { aafplot_mutations } from './modules/05-plots.nf'
include { aafplot_mutations_2 } from './modules/05-plots.nf'
include { aafplot_amplicons } from './modules/05-plots.nf'

include { ampliconsorting_DeltaReads } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_OmicronReads } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_samtools } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_bgzip } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_fasta } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_lineageAssignment_Pangolin } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_lineageAssignment_Nextclade } from './modules/06-ampliconSorting.nf'

include { report } from './modules/07-report.nf'
include { report_no_flag } from './modules/07-report.nf'


// import subworkflows
include { bammix_flagged } from './workflows/bammix_flagged.nf'
include { no_bammix_flag } from './workflows/no_bammix_flag.nf'
include { conditional } from './workflows/conditional.nf'


workflow {

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

        // Input bam, fastq, and fasta files
        ch_bam_file.map { bamfilePath -> tuple(bamfilePath) }
        ch_cat_fasta = ch_fasta.collectFile(name: 'all_sequences.fasta', newLine: true )
        ch_fastq.map { fastqPath -> tuple(fastqPath) }

        // Lineage assignment
           // through Pangolin and Nextclade using FASTA files
           pangolin( ch_cat_fasta )
           nextclade( ch_cat_fasta, params.SC2_dataset )
           lineage_assignment( pangolin.out.pangolin_csv, nextclade.out.nextclade_tsv )
           // through VirStrain using FASTQ files
           virstrain(params.virstrain_database, ch_fastq)
           virstrain_summary(params.virstrain_txt_dir, virstrain.out.txt)

        // Lineage abundance estimation per sample through Freyja using BAM files
           freyja(params.reference, ch_bam_file)
           freyja_demix(freyja.out.freyja_variants)
           freyja_aggregate(freyja_demix.out.tsv_demix.collect())
           freyja_plot(freyja_aggregate.out.freyja_aggregated_file)

        // Detecting of nucleotide mixtures from all samples
           bammix( nextclade.out.nextclade_tsv, ch_bam_file, ch_bam_index )
           bam_filter( bammix.out.bammixflagged_csv )
           // Count number of samples flagged for nucleotide mixtures
           ch_count = bam_filter.out.flatMap { it.split("\n") }.count() 

        // IF ch_count > 0 ; There are samples flagged by bammix
          // Filter high quality reads from samples with nucleotide mixture and make VCF
           ch_bammix_flagged = bam_filter.out.flatMap{ it.split("\n") }
           makevcf_2( ch_bammix_flagged, params.reference )
           bcftools(makevcf_2.out.mpileup)
        
           // Alternative Allele Fraction (AAF) plotting
           bammixplot(bcftools.out.filtered_vcf)
           aafplot_mutations_2(params.mutations_table, bcftools.out.filtered_vcf)
           aafplot_amplicons(aafplot_mutations_2.out.aafplot_mut)

        // IF ch_count = 0 ; No samples flagged by bammix 
        
        // Report generation
           report(
                 lineage_assignment.out.lineageAssign_tsv,
                 bammixplot.out.bammix_plot.flatten(),
                 freyja_plot.out.freyja_plot,
                 aafplot_mutations_2.out.aafplot_mut,
                 aafplot_amplicons.out.aafplot_amp,
                 virstrain_summary.out.tsv,
                 params.report_rmd )
}