/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run nf-core/da on the merged counts produced by nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DA } from '../../../workflows/da'

workflow RNASEQ_DA {

    take:
    ch_paramsets

    main:
    DA(ch_paramsets)

    emit:
    plot_exploratory = DA.out.plot_exploratory
    plot_volcanos    = DA.out.plot_volcanos
    shinyngs_data    = DA.out.shinyngs_data
    shinyngs_app     = DA.out.shinyngs_app_file
    report_html      = DA.out.report_html
    report_bundle    = DA.out.report_bundle
}
