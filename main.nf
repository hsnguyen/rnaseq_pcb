#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

include { getGenomeAttribute } from 'plugin/nf-core-utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = getGenomeAttribute('fasta')
params.additional_fasta = getGenomeAttribute('additional_fasta')
params.transcript_fasta = getGenomeAttribute('transcript_fasta')
params.gff              = getGenomeAttribute('gff')
params.gtf              = getGenomeAttribute('gtf')
params.gene_bed         = getGenomeAttribute('bed12')
params.bbsplit_index    = getGenomeAttribute('bbsplit')
params.sortmerna_index  = getGenomeAttribute('sortmerna')
params.star_index       = getGenomeAttribute('star')
params.rsem_index       = getGenomeAttribute('rsem')
params.hisat2_index     = getGenomeAttribute('hisat2')
params.salmon_index     = getGenomeAttribute('salmon')
params.kallisto_index   = getGenomeAttribute('kallisto')
params.bowtie2_index    = getGenomeAttribute('bowtie2')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList          } from 'plugin/nf-schema'
include { RNASEQ                     } from './workflows/rnaseq'
include { PREPARE_GENOME_REFERENCES  } from './subworkflows/local/prepare_genome_references'
include { PREPARE_GENOME_INDICES     } from './subworkflows/local/prepare_genome_indices'
include { PIPELINE_INITIALISATION    } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION        } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize         } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { defineQcTools              } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { isStarIndexLegacy          } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { RNASEQ_DA } from './subworkflows/local/da'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_RNASEQ {

    main:

    //
    // SUBWORKFLOW: Prepare reference genome files (FASTA, GTF, BED, transcript FASTA, chrom.sizes, rRNA FASTAs, Kraken DB)
    //
    PREPARE_GENOME_REFERENCES (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.ribo_database_manifest,
        params.kraken_db,
        params.gencode,
        params.gffread_transcript_fasta,
        params.featurecounts_group_type,
        params.aligner,
        params.pseudo_aligner,
        params.skip_gtf_filter,
        params.remove_ribo_rna && !(params.ribo_removal_tool == "bowtie2" && params.bowtie2_rrna_index) ? params.ribo_removal_tool : null,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.contaminant_screening,
        params.prokaryotic ?: false
    )

    //
    // SUBWORKFLOW: Build or load aligner / pseudo-aligner / filtering indices
    //
    PREPARE_GENOME_INDICES (
        PREPARE_GENOME_REFERENCES.out.fasta_fai,
        PREPARE_GENOME_REFERENCES.out.gtf,
        PREPARE_GENOME_REFERENCES.out.transcript_fasta,
        PREPARE_GENOME_REFERENCES.out.rrna_fastas,
        params.fasta ? true : false,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bowtie2_index,
        params.bbsplit_index,
        params.sortmerna_index,
        params.bowtie2_rrna_index,
        params.aligner,
        params.pseudo_aligner,
        params.skip_bbsplit,
        params.remove_ribo_rna ? params.ribo_removal_tool : null,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.use_parabricks_star,
        isStarIndexLegacy() ?: false,
        params.hisat2_build_memory,
        anySampleAutoStrandedness()
    )

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME_REFERENCES
            .out
            .fasta_fai
            .map { _meta, _fasta, fai -> checkMaxContigSize(fai) }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))
    def qc_tools = defineQcTools(params)

    RNASEQ (
        ch_samplesheet,
        PREPARE_GENOME_REFERENCES.out.fasta_fai,
        PREPARE_GENOME_REFERENCES.out.gtf,
        PREPARE_GENOME_REFERENCES.out.chrom_sizes,
        PREPARE_GENOME_REFERENCES.out.gene_bed,
        PREPARE_GENOME_REFERENCES.out.transcript_fasta,
        PREPARE_GENOME_INDICES.out.star_index,
        PREPARE_GENOME_INDICES.out.rsem_index,
        PREPARE_GENOME_INDICES.out.hisat2_index,
        PREPARE_GENOME_INDICES.out.bowtie2_index,
        PREPARE_GENOME_INDICES.out.salmon_index,
        PREPARE_GENOME_INDICES.out.kallisto_index,
        PREPARE_GENOME_INDICES.out.bbsplit_index,
        PREPARE_GENOME_REFERENCES.out.rrna_fastas,
        PREPARE_GENOME_INDICES.out.sortmerna_index,
        PREPARE_GENOME_INDICES.out.bowtie2_rrna_index,
        PREPARE_GENOME_INDICES.out.splicesites,
        PREPARE_GENOME_REFERENCES.out.kraken_db,
        qc_tools,
        params.fail_on_low_trimmed_reads,
        params.fail_on_low_mapped_reads,
        params.fail_on_undetermined_strandedness
    )

    emit:
    trim_status    = RNASEQ.out.trim_status    // channel: [id, boolean]
    map_status     = RNASEQ.out.map_status     // channel: [id, boolean]
    strand_status  = RNASEQ.out.strand_status  // channel: [id, boolean]
    multiqc_report = RNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    counts_gene_length_scaled = RNASEQ.out.counts_gene_length_scaled
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    if (params.run_da) {
        if (!params.da_contrasts || !params.da_gtf) {
            error 'When --run_da is enabled, --da_contrasts, and --da_gtf are required'
        }
        if (params.skip_quantification_merge) {
            error 'DA requires merged gene counts; do not combine --run_da with --skip_quantification_merge'
        }
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNASEQ ()

    if (params.run_da) {
        ch_da_paramsets = NFCORE_RNASEQ.out.counts_gene_length_scaled
            .map { _meta, counts ->
                [
                    id: params.da_study_name,
                    paramset_name: params.da_study_name,
                    params: [
                        study_name: params.da_study_name,
                        study_type: 'rnaseq',
                        input: file(params.da_input?:params.input, checkIfExists: true),
                        matrix: counts,
                        contrasts: file(params.da_contrasts, checkIfExists: true),
                        //round_digits: 3, 
                        outdir: "${params.outdir}/da",
                        gtf: file(params.da_gtf, checkIfExists: true),
                        gene_sets_files: (params.da_gene_sets ?: null), //can be comma-separated list, check later
                        report_file: file(params.da_report, checkIfExists: true),
                        logo_file: file(params.da_logo, checkIfExists: true),
                        css_file: file(params.da_css, checkIfExists: true),
                        citations_file: file(params.da_citations, checkIfExists: true),
                        shinyngs_build_app: params.da_shinyngs_build_app,
                        shinyngs_deploy_to_shinyapps_io: params.da_shinyngs_deploy_to_shinyapps_io,
                        shinyngs_shinyapps_account: params.da_shinyngs_shinyapps_account,
                        shinyngs_shinyapps_app_name: params.da_shinyngs_shinyapps_app_name,
                        functional_method: params.da_functional_method,
                        features_id_col: params.da_features_id_col,
                        features_name_col: params.da_features_name_col,
                        features_gtf_feature_type: params.da_features_gtf_feature_type,
                        features_gtf_table_first_field: params.da_features_gtf_table_first_field,
                        observations_id_col: params.da_observations_id_col,
                        observations_type: 'sample',
                        differential_fc_column: 'log2FoldChange',
                        differential_pval_column: 'pvalue',
                        differential_qval_column: 'padj',
                        differential_foldchanges_logged: true,
                        skip_reports: false,
                        filtering_min_samples: params.da_filtering_min_samples,
                        filtering_min_abundance: params.da_filtering_min_abundance,
                        filtering_min_proportion: params.da_filtering_min_proportion,
                        filtering_grouping_var: params.da_filtering_grouping_var,
                        filtering_min_proportion_not_na: params.da_filtering_min_proportion_not_na,
                        filtering_min_samples_not_na: params.da_filtering_min_samples_not_na,
                        exploratory_main_variable: params.da_exploratory_main_variable,
                        exploratory_clustering_method: params.da_exploratory_clustering_method,
                        exploratory_cor_method: params.da_exploratory_cor_method,
                        exploratory_n_features: params.da_exploratory_n_features,
                        exploratory_whisker_distance: params.da_exploratory_whisker_distance,
                        exploratory_mad_threshold: params.da_exploratory_mad_threshold,
                        exploratory_assay_names: params.da_exploratory_assay_names,
                        exploratory_final_assay: params.da_exploratory_final_assay,
                        exploratory_log2_assays: params.da_exploratory_log2_assays,
                        exploratory_palette_name: params.da_exploratory_palette_name,
                        differential_method: params.da_differential_method, 
                        differential_file_suffix: params.da_differential_file_suffix,
                        differential_feature_id_column: params.da_differential_feature_id_column,
                        differential_min_fold_change: params.da_differential_min_fold_change,
                        differential_max_pval: params.da_differential_max_pval,
                        differential_max_qval: params.da_differential_max_qval,
                        differential_palette_name: params.da_differential_palette_name,
                        differential_subset_to_contrast_samples: params.da_differential_subset_to_contrast_samples,
                        deseq2_test: params.da_deseq2_test,
                        deseq2_fit_type: params.da_deseq2_fit_type,
                        deseq2_sf_type: params.da_deseq2_sf_type,
                        deseq2_min_replicates_for_replace: params.da_deseq2_min_replicates_for_replace,
                        deseq2_use_t: params.da_deseq2_use_t,
                        deseq2_lfc_threshold: params.da_deseq2_lfc_threshold,
                        deseq2_alt_hypothesis: params.da_deseq2_alt_hypothesis,
                        deseq2_independent_filtering: params.da_deseq2_independent_filtering,
                        deseq2_p_adjust_method: params.da_deseq2_p_adjust_method,
                        deseq2_alpha: params.da_deseq2_alpha,
                        deseq2_minmu: params.da_deseq2_minmu,
                        deseq2_vs_method: params.da_deseq2_vs_method,
                        deseq2_shrink_lfc: params.da_deseq2_shrink_lfc,
                        deseq2_vs_blind: params.da_deseq2_vs_blind,
                        deseq2_vst_nsub: params.da_deseq2_vst_nsub,
                        gsea_nperm: params.da_gsea_nperm,
                        gsea_permute: params.da_gsea_permute,
                        gsea_scoring_scheme: params.da_gsea_scoring_scheme,
                        gsea_metric: params.da_gsea_metric,
                        gsea_sort: params.da_gsea_sort,
                        gsea_order: params.da_gsea_order,
                        gsea_set_max: params.da_gsea_set_max,
                        gsea_set_min: params.da_gsea_set_min,
                        gsea_norm: params.da_gsea_norm,
                        gsea_rnd_type: params.da_gsea_rnd_type,
                        gsea_make_sets: params.da_gsea_make_sets,
                        gsea_median: params.da_gsea_median,
                        gsea_num: params.da_gsea_num,
                        gsea_plot_top_x: params.da_gsea_plot_top_x,
                        gsea_save_rnd_lists: params.da_gsea_save_rnd_lists,
                        gprofiler2_organism: params.da_gprofiler2_organism,
                        gprofiler2_significant: params.da_gprofiler2_significant,
                        gprofiler2_measure_underrepresentation: params.da_gprofiler2_measure_underrepresentation,
                        gprofiler2_correction_method: params.da_gprofiler2_correction_method,
                        gprofiler2_sources: params.da_gprofiler2_sources,
                        gprofiler2_evcodes: params.da_gprofiler2_evcodes,
                        gprofiler2_max_qval: params.da_gprofiler2_max_qval,
                        gprofiler2_token: params.da_gprofiler2_token,
                        gprofiler2_background_file: params.da_gprofiler2_background_file,
                        gprofiler2_background_column: params.da_gprofiler2_background_column,
                        gprofiler2_domain_scope: params.da_gprofiler2_domain_scope,
                        gprofiler2_min_diff: params.da_gprofiler2_min_diff,
                        gprofiler2_palette_name: params.da_gprofiler2_palette_name,
                        decoupler_network: params.da_decoupler_network,
                        decoupler_min_n: params.da_decoupler_min_n,
                        decoupler_methods: params.da_decoupler_methods,
                        grea_set_min: params.da_grea_set_min,
                        grea_set_max: params.da_grea_set_max,
                        grea_permutation: params.da_grea_permutation,                        
                        report_title: 'PCB RNASeq Report',
                        report_author: 'Son Nguyen',
                        report_description: 'Draft report for RNASeq analysis'
                    ]
                ]
            }
        RNASEQ_DA(ch_da_paramsets)
    }

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_RNASEQ.out.multiqc_report,
        NFCORE_RNASEQ.out.trim_status,
        NFCORE_RNASEQ.out.map_status,
        NFCORE_RNASEQ.out.strand_status
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check whether any sample that will be processed from FASTQ declares strandedness 'auto'
// (strandedness is only inferred on the FASTQ path, so BAM rows under --skip_alignment do not count)
//

def anySampleAutoStrandedness() {
    samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
        .any { row -> row[0].strandedness == 'auto' && !(params.skip_alignment && (row[3] || row[4])) }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
