// ============================================================
// modules/read_retention.nf
// Aggregate per-sample read counts across all filtering stages
// ============================================================

process READ_RETENTION_SUMMARY {
    tag "${sample_id}"
    label 'process_single'

    publishDir "${params.outdir}/05_filtering/read_retention",
        mode: 'copy',
        pattern: "*.read_retention.tsv"

    input:
    tuple val(sample_id),
          path(flagstat),
          path(mito_stats),
          path(dup_metrics)   // may be null/empty if dedup was skipped

    output:
    tuple val(sample_id), path("${sample_id}.read_retention.tsv"), emit: tsv
    path "versions.yml",                                           emit: versions

    script:
    def dup_arg = dup_metrics ? "--dup_metrics \"${dup_metrics}\"" : ""
    """
    read_retention_summary.py \\
        --sample_id    "${sample_id}" \\
        --flagstat     "${flagstat}" \\
        --mito_stats   "${mito_stats}" \\
        ${dup_arg} \\
        --output       "${sample_id}.read_retention.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.read_retention.tsv versions.yml
    """
}
