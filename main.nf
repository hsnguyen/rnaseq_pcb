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
                        round_digits: 3, 
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
                        observations_id_col: params.da_observations_id_col,
                        differential_fc_column: 'log2FoldChange',
                        differential_pval_column: 'pvalue',
                        differential_qval_column: 'padj',
                        differential_foldchanges_logged: true,
                        differential_subset_to_contrast_samples: true,
                        differential_palette_name: 'Set1',
                        observations_type: 'sample',
                        differential_method: 'deseq2',
                        skip_reports: false,
                        features_gtf_feature_type: 'transcript',
                        features_gtf_table_first_field: 'gene_id',
                        filtering_min_samples: 1,
                        filtering_min_abundance: false,
                        filtering_min_proportion: 0,
                        filtering_grouping_var: null,
                        filtering_min_proportion_not_na: 0.5,
                        filtering_min_samples_not_na: 0,
                        //DESeq2 options
                        deseq2_test: 'Wald',
                        deseq2_fit_type: 'parametric',
                        deseq2_sf_type: 'ratio', 
                        deseq2_min_replicates_for_replace: 7, 
                        deseq2_use_t: false,
                        deseq2_lfc_threshold: 0, 
                        deseq2_alt_hypothesis: 'greaterAbs', 
                        deseq2_independent_filtering: true, 
                        deseq2_p_adjust_method: 'BH', 
                        deseq2_alpha: 0.1, 
                        deseq2_minmu: 0.5, 
                        deseq2_vs_method: 'vst', 
                        deseq2_vs_blind: true, 
                        deseq2_vst_nsub: 1000, 
                        deseq2_shrink_lfc: true, 
                        //GSEA
                        gsea_nperm: 1000,
                        gsea_permute: 'phenotype',
                        gsea_scoring_scheme: 'weighted',
                        gsea_metric: 'Signal2Noise',
                        gsea_sort: 'real',
                        gsea_order: 'descending',
                        gsea_set_max: 500,
                        gsea_set_min: 15,
                        gsea_norm: 'meandiv',
                        gsea_rnd_type: 'no_balance',
                        gsea_make_sets: true,
                        gsea_median: false,
                        gsea_num: 100,
                        gsea_plot_top_x: 20,
                        gsea_save_rnd_lists: false,
                        gsea_zip_report: false,
                        //GPROFILER2
                        gprofiler2_significant: true,
                        gprofiler2_measure_underrepresentation: false,
                        gprofiler2_correction_method: 'gSCS',
                        gprofiler2_evcodes: false,
                        gprofiler2_max_qval: 0.05,
                        gprofiler2_domain_scope: 'annotated',
                        gprofiler2_min_diff: 1,
                        gprofiler2_palette_name: 'Blues'
                        gprofiler2_token: null,
                        gprofiler2_organism: 'hsapiens',
                        gprofiler2_background_column: null,
                        gprofiler2_sources: null,
                        //others
                        exploratory_assay_names: 'raw,variance_stabilised',
                        exploratory_final_assay: 'variance_stabilised',
                        exploratory_mad_threshold: -5,
                        exploratory_log2_assays: 'raw,normalised'
                        reference: '',
                        target: '',
                        report_title: 'PCB RNASeq Report',
                        report_author: 'Son Nguyen',
                        report_description: 'Draft report for RNASeq analysis'

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
