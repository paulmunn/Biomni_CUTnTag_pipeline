// ============================================================
// modules/make_mqc_custom.nf
// Generate MultiQC-compatible custom content files
// ============================================================

process MAKE_MQC_CUSTOM {
    tag "make_mqc_custom"
    label 'process_single'

    publishDir "${params.outdir}/13_multiqc/custom_content",
        mode: 'copy'

    input:
    path(retention_files, stageAs: "retention/*")
    path(frip_files,      stageAs: "frip/*")
    path(peak_stats,      stageAs: "peaks/*")
    path(flagstat_files,  stageAs: "flagstat/*")

    output:
    path "read_retention_mqc.tsv",    emit: retention_mqc
    path "frip_scores_mqc.tsv",       emit: frip_mqc
    path "peak_counts_mqc.tsv",       emit: peaks_mqc
    path "mito_fraction_mqc.tsv",     emit: mito_mqc
    path "*.tsv",                     emit: all_mqc
    path "versions.yml",              emit: versions

    script:
    """
    make_mqc_custom.py \\
        --retention_dir  retention/ \\
        --frip_dir       frip/ \\
        --peaks_dir      peaks/ \\
        --flagstat_dir   flagstat/ \\
        --outdir         .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch read_retention_mqc.tsv frip_scores_mqc.tsv peak_counts_mqc.tsv mito_fraction_mqc.tsv
    touch versions.yml
    """
}
