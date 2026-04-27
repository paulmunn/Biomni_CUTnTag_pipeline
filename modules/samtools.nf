// ============================================================
// modules/samtools.nf
// SAMtools: sort, index, flagstat, idxstats, stats
// ============================================================

process SAMTOOLS_SORT_INDEX {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/04_alignment",
        mode: 'copy',
        pattern: "*.{bam,bai}"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam_bai
    path "versions.yml", emit: versions

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m 2G \\
        -o "${sample_id}_sorted.bam" \\
        "${bam}"

    samtools index \\
        -@ ${task.cpus} \\
        "${sample_id}_sorted.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_sorted.bam ${sample_id}_sorted.bam.bai versions.yml
    """
}

// ─────────────────────────────────────────────────────────────

process SAMTOOLS_FLAGSTAT {
    tag "${sample_id}"
    label 'process_single'

    publishDir "${params.outdir}/04_alignment/stats",
        mode: 'copy',
        pattern: "*.flagstat"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.flagstat"), emit: stats
    path "versions.yml",                                 emit: versions

    script:
    """
    samtools flagstat \\
        -@ ${task.cpus} \\
        "${bam}" \\
        > "${sample_id}.flagstat"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.flagstat versions.yml
    """
}

// ─────────────────────────────────────────────────────────────

process SAMTOOLS_IDXSTATS {
    tag "${sample_id}"
    label 'process_single'

    publishDir "${params.outdir}/04_alignment/stats",
        mode: 'copy',
        pattern: "*.idxstats"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.idxstats"), emit: stats
    path "versions.yml",                                  emit: versions

    script:
    """
    samtools idxstats "${bam}" > "${sample_id}.idxstats"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.idxstats versions.yml
    """
}

// ─────────────────────────────────────────────────────────────

process SAMTOOLS_STATS {
    tag "${sample_id}"
    label 'process_single'

    publishDir "${params.outdir}/04_alignment/stats",
        mode: 'copy',
        pattern: "*.stats"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.stats"), emit: stats
    path "versions.yml",                               emit: versions

    script:
    """
    samtools stats \\
        -@ ${task.cpus} \\
        "${bam}" \\
        > "${sample_id}.stats"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.stats versions.yml
    """
}
