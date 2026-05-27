// ============================================================
// modules/filter.nf
// Filtering: mitochondrial removal, blacklist, deduplication
// ============================================================

// ── 1. Remove mitochondrial reads ────────────────────────────
process FILTER_MITO {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/05_filtering/mito_removed",
        mode: 'copy',
        pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/05_filtering/stats",
        mode: 'copy',
        pattern: "*.mito_stats.txt"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_nomito.bam"), path("${sample_id}_nomito.bam.bai"), emit: bam_bai
    tuple val(sample_id), path("${sample_id}.mito_stats.txt"),                                  emit: stats
    path "versions.yml",                                                                         emit: versions

    script:
    def mito_pattern = params.mito_name ?: "chrM|MT|M|mitochondria|Mito"
    """
    # Count total reads before mito removal
    TOTAL_READS=\$(samtools view -c -F 4 "${bam}")

    # Count mitochondrial reads
    MITO_READS=\$(samtools view -c "${bam}" \$(samtools view -H "${bam}" | \\
        grep "^@SQ" | awk '{print \$2}' | sed 's/SN://' | \\
        grep -E "^(${mito_pattern})\$" | tr '\\n' ' ') 2>/dev/null || echo 0)

    # Get list of non-mitochondrial chromosomes
    NON_MITO_CHROMS=\$(samtools view -H "${bam}" | \\
        grep "^@SQ" | awk '{print \$2}' | sed 's/SN://' | \\
        grep -vE "^(${mito_pattern})\$" | tr '\\n' ' ')

    # Filter: keep only non-mito, properly paired, MAPQ >= threshold
    samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -f 2 \\
        -q ${params.mapq} \\
        -o "${sample_id}_nomito.bam" \\
        "${bam}" \\
        \${NON_MITO_CHROMS}

    samtools index -@ ${task.cpus} "${sample_id}_nomito.bam"

    # Count reads after mito removal
    NOMITO_READS=\$(samtools view -c -F 4 "${sample_id}_nomito.bam")

    # Write stats
    echo -e "sample_id\\ttotal_reads\\tmito_reads\\tnomito_reads\\tmito_fraction" \\
        > "${sample_id}.mito_stats.txt"
    echo -e "${sample_id}\\t\${TOTAL_READS}\\t\${MITO_READS}\\t\${NOMITO_READS}\\t\$(echo "scale=4; \${MITO_READS}/\${TOTAL_READS}" | bc)" \\
        >> "${sample_id}.mito_stats.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_nomito.bam ${sample_id}_nomito.bam.bai
    echo -e "sample_id\\ttotal_reads\\tmito_reads\\tnomito_reads\\tmito_fraction" > ${sample_id}.mito_stats.txt
    echo -e "${sample_id}\\t1000000\\t50000\\t950000\\t0.05" >> ${sample_id}.mito_stats.txt
    touch versions.yml
    """
}

// ── 2. Remove blacklisted regions ────────────────────────────
process FILTER_BLACKLIST {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/05_filtering/blacklist_removed",
        mode: 'copy',
        pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/05_filtering/stats",
        mode: 'copy',
        pattern: "*.blacklist_stats.txt"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path  blacklist_bed

    output:
    tuple val(sample_id), path("${sample_id}_noblacklist.bam"), path("${sample_id}_noblacklist.bam.bai"), emit: bam_bai
    tuple val(sample_id), path("${sample_id}.blacklist_stats.txt"),                                        emit: stats
    path "versions.yml",                                                                                   emit: versions

    script:
    """
    BEFORE=\$(samtools view -c -F 4 "${bam}")

    bedtools intersect \\
        -v \\
        -abam "${bam}" \\
        -b "${blacklist_bed}" \\
    | samtools sort -@ ${task.cpus} -o "${sample_id}_noblacklist.bam"

    samtools index -@ ${task.cpus} "${sample_id}_noblacklist.bam"

    AFTER=\$(samtools view -c -F 4 "${sample_id}_noblacklist.bam")
    REMOVED=\$(( BEFORE - AFTER ))

    echo -e "sample_id\\tbefore_blacklist\\tafter_blacklist\\tremoved_blacklist" \\
        > "${sample_id}.blacklist_stats.txt"
    echo -e "${sample_id}\\t\${BEFORE}\\t\${AFTER}\\t\${REMOVED}" \\
        >> "${sample_id}.blacklist_stats.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_noblacklist.bam ${sample_id}_noblacklist.bam.bai
    touch ${sample_id}.blacklist_stats.txt versions.yml
    """
}

// ── 3a. Picard MarkDuplicates ─────────────────────────────────
// Recommendation: Picard MarkDuplicates is the standard for
// CUT&Tag. It handles paired-end optical and PCR duplicates.
// For UMI-based libraries, use UMI_TOOLS_DEDUP instead.

process PICARD_MARKDUPLICATES {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/05_filtering/dedup",
        mode: 'copy',
        pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/05_filtering/stats",
        mode: 'copy',
        pattern: "*.dup_metrics.txt"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_dedup.bam"), path("${sample_id}_dedup.bam.bai"), emit: bam_bai
    tuple val(sample_id), path("${sample_id}.dup_metrics.txt"),                               emit: metrics
    path "versions.yml",                                                                       emit: versions

    script:
    def avail_mem = (task.memory.toGiga() * 0.8).intValue()
    """
    # picard -Xmx${avail_mem}g MarkDuplicates \\
    java -jar /programs/picard-tools-3.4.0/picard.jar MarkDuplicates \\
        INPUT="${bam}" \\
        OUTPUT="${sample_id}_dedup.bam" \\
        METRICS_FILE="${sample_id}.dup_metrics.txt" \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORTED=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=./tmp_picard

    samtools index -@ ${task.cpus} "${sample_id}_dedup.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -i version | head -1 | sed 's/.*Version://')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_dedup.bam ${sample_id}_dedup.bam.bai ${sample_id}.dup_metrics.txt versions.yml
    """
}

// ── 3b. UMI-tools deduplication (optional) ───────────────────
process UMI_TOOLS_DEDUP {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/05_filtering/dedup",
        mode: 'copy',
        pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/05_filtering/stats",
        mode: 'copy',
        pattern: "*.umi_dedup.log"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_dedup.bam"), path("${sample_id}_dedup.bam.bai"), emit: bam_bai
    tuple val(sample_id), path("${sample_id}.umi_dedup.log"),                                 emit: metrics
    path "versions.yml",                                                                       emit: versions

    script:
    """
    umi_tools dedup \\
        --paired \\
        --stdin="${bam}" \\
        --stdout="${sample_id}_dedup.bam" \\
        --log="${sample_id}.umi_dedup.log" \\
        --umi-separator="${params.umi_separator}"

    samtools index -@ ${task.cpus} "${sample_id}_dedup.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: \$(umi_tools --version | sed 's/UMI-tools version: //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_dedup.bam ${sample_id}_dedup.bam.bai ${sample_id}.umi_dedup.log versions.yml
    """
}
