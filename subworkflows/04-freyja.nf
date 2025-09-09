// enable dsl2
nextflow.enable.dsl=2

// import modules
include { variants } from '../modules/04-freyja.nf'
include { demix } from '../modules/04-freyja.nf'
include { aggregate } from '../modules/04-freyja.nf'
include { plot_summarized } from '../modules/04-freyja.nf'
include { plot_lineage } from '../modules/04-freyja.nf'
include { list_lineages } from '../modules/04-freyja.nf'
include { get_lineage_def } from '../modules/04-freyja.nf'
include { mutations } from '../modules/04-freyja.nf'

workflow freyja {
    take:
        ch_bam_file
        flagged

    main:
        variants(params.reference, ch_bam_file)
        demix(variants.out.freyja_variants)
        aggregate(demix.out.tsv_demix.collect())
        plot_summarized(aggregate.out.freyja_aggregated_file)
        plot_lineage(plot_summarized.out.aggregated_tsv)

    // Only for samples flagged by bammix
        ch_freyja_tsv = demix.out.tsv_demix.collect()
        ch_freyja_samples = ch_freyja_tsv.flatten()
            .map { tsv ->
               def sample_name = tsv.baseName
               tuple(sample_name, tsv)
            }


        ch_freyja_flagged_tsv = flagged
            .join (ch_freyja_samples)
            .map { sample, tsv -> tuple (sample, tsv) }
                  
        list_lineages(ch_freyja_flagged_tsv)
        get_lineage_def(list_lineages.out.freyja_list_lin, params.annot, params.ref)
        freyja_result = get_lineage_def.out.lin_mut_tsv

        mutations(freyja_result)

    emit:
        mutations_tsv = mutations.out.processed_mut_tsv
        freyja_summarized = plot_summarized.out.png
        freyja_lineage = plot_lineage.out.png
        freyja_result = get_lineage_def.out.lin_mut_tsv
}