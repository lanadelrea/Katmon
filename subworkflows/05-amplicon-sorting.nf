// enable dsl2
nextflow.enable.dsl=2

// import modules
include { get_pos_mut } from '../modules/05-ampliconSorting.nf'
include { ampliconsorting } from '../modules/05-ampliconSorting.nf'
include { consensus } from '../modules/05-ampliconSorting.nf'
include { renamefasta } from '../modules/05-ampliconSorting.nf'
include { pangolin } from '../modules/05-ampliconSorting.nf'
include { nextclade } from '../modules/05-ampliconSorting.nf'
include { summary } from '../modules/05-ampliconSorting.nf'

workflow amplicon_sorting {
   take:
       filtered_vcf
       freyja_result

   main:
       // Amplicon sorting
       filtered_vcf
          .join(freyja_result) // only get mutations list for bammix flagged samples
          .set { get_posmut_flagged }
       get_pos_mut( get_posmut_flagged ) // List mutations for the two most abundant lineages

       get_pos_mut.out.pos_mut_lineage_A
          .join(get_pos_mut.out.pos_mut_lineage_B)
          .set { pos_mut }
       ampliconsorting(pos_mut, params.jvarkit_jar, params.sort_reads)

       ampliconsorting.out.sorted_bam_lineage_A
          .join(ampliconsorting.out.sorted_bam_lineage_B)
          .set { ampsort_bam }
       consensus(ampsort_bam, params.reference)

       consensus.out.consensus_lineage_A
          .join(consensus.out.consensus_lineage_B)
          .set { ampsort_consensus }
       renamefasta(ampsort_consensus)
       renamefasta.out.fasta
          .collectFile(name: 'all_sequences_ampsort.fasta', newLine: true )
          .set { ch_ampsort_cat }

       pangolin(ch_ampsort_cat)
       nextclade(ch_ampsort_cat, params.SC2_dataset)

       summary(pangolin.out.csv, nextclade.out.tsv)
   
   emit:
       table = summary.out.table


}