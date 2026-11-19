/*
========================================================================================
    Run DESeq2 and filter differential results for RNA-seq counts
========================================================================================
*/

include { DESEQ2_DIFFERENTIAL } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM } from '../../../modules/nf-core/deseq2/differential/main'
include { CUSTOM_FILTERDIFFERENTIALTABLE } from '../../../modules/nf-core/custom/filterdifferentialtable/main'

workflow ABUNDANCE_DIFFERENTIAL_FILTER {
    take:
    ch_input
    ch_samplesheet
    ch_transcript_lengths
    ch_control_features
    ch_contrasts

    main:
    def criteria = multiMapCriteria { meta, abundance, _method, fc_threshold, stat_threshold, samplesheet, transcript_length, control_features, meta_contrast, variable, reference, target, formula, comparison ->
        def meta_with_method = meta + [differential_method: 'deseq2']
        def meta_for_diff = meta_contrast + meta_with_method
        samples_and_matrix:
            [meta_with_method, samplesheet, abundance]
        contrasts_for_diff:
            [meta_for_diff, variable, reference, target, formula, comparison]
        filter_params:
            [meta_for_diff, [fc_threshold: fc_threshold, stat_threshold: stat_threshold]]
        transcript_length:
            [meta_with_method, (transcript_length) ?: []]
        control_features:
            [meta_with_method, (control_features) ?: []]
    }

    ch_samplesheet_with_control = ch_samplesheet
        .join(ch_transcript_lengths, remainder: true)
        .join(ch_control_features, remainder: true)

    inputs = ch_input
        .combine(ch_samplesheet_with_control, by: 0)
        .combine(ch_contrasts.transpose(), by: 0)
        .multiMap(criteria)

    norm_inputs = ch_input
        .combine(ch_samplesheet_with_control, by: 0)
        .combine(ch_contrasts.map { meta, contrasts, variables, references, targets, formulas, comparisons ->
            tuple(meta, contrasts[0], variables[0], references[0], targets[0], formulas[0], comparisons[0])
        }, by: 0)
        .multiMap(criteria)

    DESEQ2_NORM(
        norm_inputs.contrasts_for_diff,
        norm_inputs.samples_and_matrix,
        norm_inputs.control_features,
        norm_inputs.transcript_length
    )

    DESEQ2_DIFFERENTIAL(
        inputs.contrasts_for_diff,
        inputs.samples_and_matrix,
        inputs.control_features,
        inputs.transcript_length
    )

    ch_diff_filter_params = DESEQ2_DIFFERENTIAL.out.results
        .join(inputs.filter_params)
        .multiMap { meta, results, filter_meta ->
            filter_input: [meta + filter_meta, results]
            fc_input: ['log2FoldChange', filter_meta.fc_threshold, '>=']
            stat_input: ['padj', filter_meta.stat_threshold, '<=']
        }

    CUSTOM_FILTERDIFFERENTIALTABLE(
        ch_diff_filter_params.filter_input,
        ch_diff_filter_params.fc_input,
        ch_diff_filter_params.stat_input
    )

    emit:
    results_genewise           = DESEQ2_DIFFERENTIAL.out.results
    results_genewise_filtered  = CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered
    normalised_matrix          = DESEQ2_NORM.out.normalised_counts
    variance_stabilised_matrix = DESEQ2_NORM.out.rlog_counts.mix(DESEQ2_NORM.out.vst_counts).groupTuple()
    model                      = DESEQ2_DIFFERENTIAL.out.model
    size_factors               = DESEQ2_DIFFERENTIAL.out.size_factors
    dispersion_plot            = DESEQ2_DIFFERENTIAL.out.dispersion_plot_png
    md_plot                    = channel.empty()
    rdata                      = DESEQ2_DIFFERENTIAL.out.rdata
    session_info               = DESEQ2_DIFFERENTIAL.out.session_info
    versions                   = DESEQ2_DIFFERENTIAL.out.versions.mix(DESEQ2_NORM.out.versions)
}
