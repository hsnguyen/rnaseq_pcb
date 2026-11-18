/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run nf-core/differentialabundance on the merged counts produced by nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DIFFERENTIALABUNDANCE } from '../../../workflows/differentialabundance'

workflow RNASEQ_DIFFERENTIALABUNDANCE {

    take:
    ch_paramsets

    main:
    DIFFERENTIALABUNDANCE(ch_paramsets)

    emit:
    plot_exploratory = DIFFERENTIALABUNDANCE.out.plot_exploratory
    plot_volcanos    = DIFFERENTIALABUNDANCE.out.plot_volcanos
    shinyngs_data    = DIFFERENTIALABUNDANCE.out.shinyngs_data
    shinyngs_app     = DIFFERENTIALABUNDANCE.out.shinyngs_app_file
    report_html      = DIFFERENTIALABUNDANCE.out.report_html
    report_bundle    = DIFFERENTIALABUNDANCE.out.report_bundle
}
