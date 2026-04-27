// ============================================================
// modules/annotate_peaks.nf
// Peak annotation using ChIPseeker (R/Bioconductor) or HOMER
// ============================================================
// Recommendation: ChIPseeker is preferred for its rich R-based
// output, ggplot2 visualizations, and support for any GTF/TxDb.
// HOMER annotatePeaks.pl is a fast alternative for large peak sets.
// Both are supported; select with --annotation_tool chipseeker|homer

process ANNOTATE_PEAKS {
    tag "${peak_id}"
    label 'process_medium'

    publishDir "${params.outdir}/10_peak_annotation/${peak_id}",
        mode: 'copy'

    input:
    tuple val(peak_id), path(peak_file)
    path  annotation_gtf

    output:
    tuple val(peak_id), path("${peak_id}_annotated.tsv"),          emit: annotated
    tuple val(peak_id), path("${peak_id}_annotation_summary.tsv"), emit: summary
    tuple val(peak_id), path("${peak_id}_annotation_pie.png"),     emit: pie_plot,  optional: true
    tuple val(peak_id), path("${peak_id}_annotation_bar.png"),     emit: bar_plot,  optional: true
    path "versions.yml",                                           emit: versions

    script:
    if (params.annotation_tool == 'homer') {
        """
        # HOMER annotation
        annotatePeaks.pl \\
            "${peak_file}" \\
            none \\
            -gtf "${annotation_gtf}" \\
            -cpu ${task.cpus} \\
            > "${peak_id}_annotated.tsv"

        # Generate summary from HOMER output
        awk 'NR>1 {print \$8}' "${peak_id}_annotated.tsv" | \\
            sort | uniq -c | sort -rn | \\
            awk 'BEGIN{print "annotation\\tcount"} {print \$2"\\t"\$1}' \\
            > "${peak_id}_annotation_summary.tsv"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            homer: \$(annotatePeaks.pl 2>&1 | head -3 | grep -i version | sed 's/.*version //')
        END_VERSIONS
        """
    } else {
        // Default: ChIPseeker
        """
        annotate_peaks_chipseeker.R \\
            --peaks    "${peak_file}" \\
            --gtf      "${annotation_gtf}" \\
            --prefix   "${peak_id}" \\
            --outdir   .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
            bioconductor-chipseeker: \$(Rscript -e "cat(as.character(packageVersion('ChIPseeker')))")
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${peak_id}_annotated.tsv ${peak_id}_annotation_summary.tsv
    touch ${peak_id}_annotation_pie.png ${peak_id}_annotation_bar.png
    touch versions.yml
    """
}
