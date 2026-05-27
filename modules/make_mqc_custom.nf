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
    path retention_files   // plain path — no stageAs; collected list staged flat
    path frip_files
    path peak_stats
    path flagstat_files

    output:
    path "read_retention_mqc.tsv",    emit: retention_mqc
    path "frip_scores_mqc.tsv",       emit: frip_mqc
    path "peak_counts_mqc.tsv",       emit: peaks_mqc
    path "mito_fraction_mqc.tsv",     emit: mito_mqc
    path "versions.yml",              emit: versions

    script:
    // Files are staged flat into the work dir; pass '.' as each dir
    // so the Python script finds them via glob in the working directory.
    """
    mkdir -p retention frip peaks flagstat

    # Move staged files into named subdirs so the Python script can find them
    for f in *.read_retention.tsv; do [ -f "\$f" ] && mv "\$f" retention/ || true; done
    for f in *.frip.txt;           do [ -f "\$f" ] && mv "\$f" frip/      || true; done
    for f in *_peak_stats.txt;     do [ -f "\$f" ] && mv "\$f" peaks/     || true; done
    for f in *.flagstat;           do [ -f "\$f" ] && mv "\$f" flagstat/  || true; done

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
