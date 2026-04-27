// ============================================================
// modules/macs3.nf
// MACS3 peak calling: per-sample and group-level
// ============================================================
// Recommendation: For CUT&Tag, use --format BAMPE (paired-end BAM)
// which uses actual fragment sizes rather than estimated extension.
// Narrow peaks (--call-summits) are appropriate for TF and sharp
// histone marks (H3K4me3, H3K27ac). Broad peaks are appropriate
// for H3K27me3, H3K9me3, H3K36me3.
// The --nomodel flag is NOT recommended for CUT&Tag; MACS3 can
// model the fragment size distribution from paired-end data.

process MACS3_CALLPEAK_SAMPLE {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/08_peaks_per_sample/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id),
          path(bam), path(bai),
          path(ctrl_bam), path(ctrl_bai),
          val(assoc)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), emit: peaks,    optional: true
    tuple val(sample_id), path("${sample_id}_peaks.broadPeak"),  emit: broad,    optional: true
    tuple val(sample_id), path("${sample_id}_summits.bed"),      emit: summits,  optional: true
    tuple val(sample_id), path("${sample_id}_peaks.xls"),        emit: xls
    path "${sample_id}_macs3.log",                               emit: log
    tuple val(sample_id), path("${sample_id}_peak_stats.txt"),   emit: stats
    path "versions.yml",                                         emit: versions

    script:
    def peak_mode  = assoc.peak_mode ?: params.macs3_genome
    def broad_flag = (peak_mode == 'broad') ? "--broad --broad-cutoff ${params.macs3_qvalue}" : "--call-summits"
    def ctrl_arg   = ctrl_bam ? "--control ${ctrl_bam}" : ""
    def qval       = params.macs3_qvalue ?: 0.05
    def genome     = params.macs3_genome ?: "hs"
    """
    macs3 callpeak \\
        --treatment "${bam}" \\
        ${ctrl_arg} \\
        --format BAMPE \\
        --gsize ${genome} \\
        --name "${sample_id}" \\
        --outdir . \\
        --qvalue ${qval} \\
        ${broad_flag} \\
        --keep-dup all \\
        --verbose 3 \\
        2> "${sample_id}_macs3.log"

    # Generate peak stats
    PEAK_FILE="${sample_id}_peaks.narrowPeak"
    if [ ! -f "\${PEAK_FILE}" ]; then
        PEAK_FILE="${sample_id}_peaks.broadPeak"
    fi

    if [ -f "\${PEAK_FILE}" ]; then
        N_PEAKS=\$(wc -l < "\${PEAK_FILE}")
        MEDIAN_SCORE=\$(awk '{print \$5}' "\${PEAK_FILE}" | sort -n | awk 'BEGIN{c=0;s=0} {a[c++]=\$1; s+=\$1} END{print (c%2==0)?(a[int(c/2)-1]+a[int(c/2)])/2:a[int(c/2)]}')
    else
        N_PEAKS=0
        MEDIAN_SCORE=0
    fi

    echo -e "sample_id\\tn_peaks\\tmedian_score\\tpeak_mode" > "${sample_id}_peak_stats.txt"
    echo -e "${sample_id}\\t\${N_PEAKS}\\t\${MEDIAN_SCORE}\\t${peak_mode}" >> "${sample_id}_peak_stats.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_peaks.narrowPeak ${sample_id}_summits.bed ${sample_id}_peaks.xls
    touch ${sample_id}_macs3.log ${sample_id}_peak_stats.txt versions.yml
    """
}

// ── Group-level peak calling ──────────────────────────────────

process MACS3_CALLPEAK_GROUP {
    tag "${merge_group_id}"
    label 'process_medium'

    publishDir "${params.outdir}/09_peaks_merged_groups/${merge_group_id}",
        mode: 'copy'

    input:
    tuple val(merge_group_id),
          path(merged_bam), path(merged_bai),
          path(ctrl_bam),   path(ctrl_bai),
          val(peak_mode)

    output:
    tuple val(merge_group_id), path("${merge_group_id}_peaks.narrowPeak"), emit: peaks,   optional: true
    tuple val(merge_group_id), path("${merge_group_id}_peaks.broadPeak"),  emit: broad,   optional: true
    tuple val(merge_group_id), path("${merge_group_id}_summits.bed"),      emit: summits, optional: true
    tuple val(merge_group_id), path("${merge_group_id}_peaks.xls"),        emit: xls
    path "${merge_group_id}_macs3.log",                                    emit: log
    tuple val(merge_group_id), path("${merge_group_id}_peak_stats.txt"),   emit: stats
    path "versions.yml",                                                   emit: versions

    script:
    def broad_flag = (peak_mode == 'broad') ? "--broad --broad-cutoff ${params.macs3_qvalue}" : "--call-summits"
    def ctrl_arg   = ctrl_bam ? "--control ${ctrl_bam}" : ""
    def qval       = params.macs3_qvalue ?: 0.05
    def genome     = params.macs3_genome ?: "hs"
    """
    macs3 callpeak \\
        --treatment "${merged_bam}" \\
        ${ctrl_arg} \\
        --format BAMPE \\
        --gsize ${genome} \\
        --name "${merge_group_id}" \\
        --outdir . \\
        --qvalue ${qval} \\
        ${broad_flag} \\
        --keep-dup all \\
        --verbose 3 \\
        2> "${merge_group_id}_macs3.log"

    PEAK_FILE="${merge_group_id}_peaks.narrowPeak"
    if [ ! -f "\${PEAK_FILE}" ]; then
        PEAK_FILE="${merge_group_id}_peaks.broadPeak"
    fi

    if [ -f "\${PEAK_FILE}" ]; then
        N_PEAKS=\$(wc -l < "\${PEAK_FILE}")
        MEDIAN_SCORE=\$(awk '{print \$5}' "\${PEAK_FILE}" | sort -n | awk 'BEGIN{c=0} {a[c++]=\$1} END{print (c%2==0)?(a[int(c/2)-1]+a[int(c/2)])/2:a[int(c/2)]}')
    else
        N_PEAKS=0
        MEDIAN_SCORE=0
    fi

    echo -e "merge_group_id\\tn_peaks\\tmedian_score\\tpeak_mode" > "${merge_group_id}_peak_stats.txt"
    echo -e "${merge_group_id}\\t\${N_PEAKS}\\t\${MEDIAN_SCORE}\\t${peak_mode}" >> "${merge_group_id}_peak_stats.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //')
    END_VERSIONS
    """

    stub:
    """
    touch ${merge_group_id}_peaks.narrowPeak ${merge_group_id}_summits.bed ${merge_group_id}_peaks.xls
    touch ${merge_group_id}_macs3.log ${merge_group_id}_peak_stats.txt versions.yml
    """
}
