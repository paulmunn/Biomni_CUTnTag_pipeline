// ============================================================
// modules/multiqc.nf
// Aggregate QC reports with MultiQC
// ============================================================

process MULTIQC {
    tag "multiqc_${stage}"
    label 'process_low'

    publishDir "${outdir_path}", mode: 'copy'

    input:
    path(qc_files, stageAs: "qc_inputs/*")
    val  stage
    val  outdir_path

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",               emit: data
    path "*_plots",              emit: plots, optional: true
    path "versions.yml",         emit: versions

    script:
    def config_arg = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    def title_arg  = params.multiqc_title  ? "--title \"${params.multiqc_title} - ${stage}\"" : ""
    """
    multiqc \\
        --force \\
        --filename "${stage}_multiqc_report.html" \\
        ${config_arg} \\
        ${title_arg} \\
        qc_inputs/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${stage}_multiqc_report.html
    mkdir -p ${stage}_multiqc_report_data
    touch versions.yml
    """
}
