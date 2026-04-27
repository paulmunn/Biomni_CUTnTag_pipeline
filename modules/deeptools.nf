// ============================================================
// modules/deeptools.nf
// deepTools: BigWig, BedGraph, computeMatrix, plotProfile, plotHeatmap
// ============================================================

// ── BigWig signal tracks ──────────────────────────────────────
// Recommendation: CPM normalization is appropriate for most
// CUT&Tag experiments. Use RPGC if comparing across experiments
// with different sequencing depths and you have a reliable
// effective genome size. Bin size of 10 bp is a good default.

process BAMCOVERAGE_BIGWIG {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/06_bigwig",
        mode: 'copy',
        pattern: "*.bigwig"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path  chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.bigwig"), emit: bigwig
    path "versions.yml",                               emit: versions

    script:
    def norm_method  = params.normalization ?: "CPM"
    def bin_size     = params.bigwig_binsize ?: 10
    def eff_gs_arg   = (norm_method == "RPGC") ? "--effectiveGenomeSize ${params.effective_genome_size}" : ""
    def blacklist_arg = (params.blacklist && !params.skip_blacklist) ? "--blackListFileName ${params.blacklist}" : ""
    """
    bamCoverage \\
        --bam "${bam}" \\
        --outFileName "${sample_id}.bigwig" \\
        --outFileFormat bigwig \\
        --binSize ${bin_size} \\
        --normalizeUsing ${norm_method} \\
        ${eff_gs_arg} \\
        ${blacklist_arg} \\
        --extendReads \\
        --ignoreDuplicates \\
        --minMappingQuality ${params.mapq} \\
        --numberOfProcessors ${task.cpus} \\
        --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed 's/bamCoverage //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.bigwig versions.yml
    """
}

// ── BedGraph signal tracks ────────────────────────────────────

process BAMCOVERAGE_BEDGRAPH {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/07_bedgraph",
        mode: 'copy',
        pattern: "*.bedgraph"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path  chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.bedgraph"), emit: bedgraph
    path "versions.yml",                                 emit: versions

    script:
    def norm_method  = params.normalization ?: "CPM"
    def bin_size     = params.bigwig_binsize ?: 10
    def eff_gs_arg   = (norm_method == "RPGC") ? "--effectiveGenomeSize ${params.effective_genome_size}" : ""
    def blacklist_arg = (params.blacklist && !params.skip_blacklist) ? "--blackListFileName ${params.blacklist}" : ""
    """
    bamCoverage \\
        --bam "${bam}" \\
        --outFileName "${sample_id}.bedgraph" \\
        --outFileFormat bedgraph \\
        --binSize ${bin_size} \\
        --normalizeUsing ${norm_method} \\
        ${eff_gs_arg} \\
        ${blacklist_arg} \\
        --extendReads \\
        --ignoreDuplicates \\
        --minMappingQuality ${params.mapq} \\
        --numberOfProcessors ${task.cpus}

    # Sort bedgraph for genome browser compatibility
    sort -k1,1 -k2,2n "${sample_id}.bedgraph" > "${sample_id}_sorted.bedgraph"
    mv "${sample_id}_sorted.bedgraph" "${sample_id}.bedgraph"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed 's/bamCoverage //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.bedgraph versions.yml
    """
}

// ── computeMatrix around TSS ──────────────────────────────────

process COMPUTE_MATRIX_TSS {
    tag "computeMatrix_tss"
    label 'process_high'

    publishDir "${params.outdir}/12_deeptools/matrices",
        mode: 'copy',
        pattern: "*.gz"
    publishDir "${params.outdir}/12_deeptools/matrices",
        mode: 'copy',
        pattern: "*.tab"

    input:
    path bigwig_files
    path tss_bed
    val  tss_window

    output:
    path "matrix_tss.gz",  emit: matrix
    path "matrix_tss.tab", emit: table, optional: true
    path "versions.yml",   emit: versions

    script:
    def bw_list = bigwig_files instanceof List ? bigwig_files.join(' ') : bigwig_files
    """
    computeMatrix reference-point \\
        --scoreFileName ${bw_list} \\
        --regionsFileName "${tss_bed}" \\
        --referencePoint TSS \\
        --beforeRegionStartLength ${tss_window} \\
        --afterRegionStartLength  ${tss_window} \\
        --binSize 10 \\
        --missingDataAsZero \\
        --skipZeros \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName matrix_tss.gz \\
        --outFileNameMatrix matrix_tss.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed 's/computeMatrix //')
    END_VERSIONS
    """

    stub:
    """
    touch matrix_tss.gz matrix_tss.tab versions.yml
    """
}

// ── computeMatrix around peak centers ────────────────────────

process COMPUTE_MATRIX_PEAKS {
    tag "computeMatrix_peaks"
    label 'process_high'

    publishDir "${params.outdir}/12_deeptools/matrices",
        mode: 'copy',
        pattern: "*.gz"

    input:
    path bigwig_files
    path peak_files
    val  window

    output:
    path "matrix_peaks.gz", emit: matrix
    path "versions.yml",    emit: versions

    script:
    def bw_list   = bigwig_files instanceof List ? bigwig_files.join(' ') : bigwig_files
    def peak_list = peak_files instanceof List ? peak_files.join(' ') : peak_files
    """
    computeMatrix reference-point \\
        --scoreFileName ${bw_list} \\
        --regionsFileName ${peak_list} \\
        --referencePoint center \\
        --beforeRegionStartLength ${window} \\
        --afterRegionStartLength  ${window} \\
        --binSize 10 \\
        --missingDataAsZero \\
        --skipZeros \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName matrix_peaks.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed 's/computeMatrix //')
    END_VERSIONS
    """

    stub:
    """
    touch matrix_peaks.gz versions.yml
    """
}

// ── plotProfile ───────────────────────────────────────────────

process PLOT_PROFILE {
    tag "plotProfile_${label}"
    label 'process_low'

    publishDir "${params.outdir}/12_deeptools/profiles",
        mode: 'copy',
        pattern: "*.{png,pdf}"

    input:
    path matrix
    val  label

    output:
    path "profile_${label}.png", emit: png
    path "profile_${label}.pdf", emit: pdf, optional: true
    path "versions.yml",         emit: versions

    script:
    """
    plotProfile \\
        --matrixFile "${matrix}" \\
        --outFileName "profile_${label}.png" \\
        --outFileSortedRegions "profile_${label}_sorted_regions.bed" \\
        --plotType lines \\
        --perGroup \\
        --dpi 300 \\
        --plotTitle "Signal Profile - ${label}"

    plotProfile \\
        --matrixFile "${matrix}" \\
        --outFileName "profile_${label}.pdf" \\
        --plotType lines \\
        --perGroup \\
        --plotTitle "Signal Profile - ${label}" || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotProfile --version | sed 's/plotProfile //')
    END_VERSIONS
    """

    stub:
    """
    touch profile_${label}.png profile_${label}.pdf versions.yml
    """
}

// ── plotHeatmap ───────────────────────────────────────────────

process PLOT_HEATMAP {
    tag "plotHeatmap_${label}"
    label 'process_low'

    publishDir "${params.outdir}/12_deeptools/heatmaps",
        mode: 'copy',
        pattern: "*.{png,pdf}"

    input:
    path matrix
    val  label

    output:
    path "heatmap_${label}.png", emit: png
    path "heatmap_${label}.pdf", emit: pdf, optional: true
    path "versions.yml",         emit: versions

    script:
    """
    plotHeatmap \\
        --matrixFile "${matrix}" \\
        --outFileName "heatmap_${label}.png" \\
        --colorMap RdBu_r \\
        --whatToShow 'heatmap and colorbar' \\
        --dpi 300 \\
        --plotTitle "Signal Heatmap - ${label}"

    plotHeatmap \\
        --matrixFile "${matrix}" \\
        --outFileName "heatmap_${label}.pdf" \\
        --colorMap RdBu_r \\
        --whatToShow 'heatmap and colorbar' \\
        --plotTitle "Signal Heatmap - ${label}" || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotHeatmap --version | sed 's/plotHeatmap //')
    END_VERSIONS
    """

    stub:
    """
    touch heatmap_${label}.png heatmap_${label}.pdf versions.yml
    """
}
