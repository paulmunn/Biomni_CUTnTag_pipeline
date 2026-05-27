// ============================================================
// modules/multiqc.nf
// Aggregate QC reports with MultiQC
// ============================================================
// Input design note:
//   qc_files is declared as `path qc_files` (no stageAs).
//   main.nf passes a .collect()-ed list of bare path values as a
//   single channel item.  Nextflow stages every file in that list
//   into the process work directory, and MultiQC is pointed at '.'
//   so it finds them all regardless of their staged names.
//   This is the standard nf-core MultiQC pattern and avoids the
//   DataflowStreamReadAdapter error caused by stageAs on a list.
//
// Output design note:
//   Outputs use explicit file names (no globs) to avoid emitting
//   DataflowStream objects.  The report and data dir names are
//   constructed from the `stage` val input.
// ============================================================

process MULTIQC {
    tag "multiqc_${stage}"
    label 'process_low'

    publishDir "${outdir_path}", mode: 'copy'

    input:
    path qc_files   // collected list of bare path values staged into work dir
    val  stage
    val  outdir_path

    output:
    path "${stage}_multiqc_report.html", emit: report
    path "${stage}_multiqc_report_data", emit: data
    path "versions.yml",                 emit: versions

    script:
    def config_arg = params.multiqc_config ? "--config '${params.multiqc_config}'" : ""
    def title_arg  = params.multiqc_title  ? "--title '${params.multiqc_title} - ${stage}'" : ""
    """
    multiqc \\
        --force \\
        --filename "${stage}_multiqc_report.html" \\
        ${config_arg} \\
        ${title_arg} \\
        .

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
