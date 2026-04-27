// ============================================================
// modules/bowtie2.nf
// Align paired-end reads with Bowtie2
// ============================================================
// Recommendation: --local --very-sensitive-local is preferred for
// CUT&Tag because reads may be slightly shorter after trimming and
// local alignment tolerates soft-clipping at read ends. The fragment
// size range -I 10 -X 700 captures nucleosome-free and mono-
// nucleosomal fragments typical of CUT&Tag.

process BOWTIE2_ALIGN {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/04_alignment/logs",
        mode: 'copy',
        pattern: "*.log"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path  bowtie2_index_dir

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}_bowtie2.log",                emit: log
    path "versions.yml",                            emit: versions

    script:
    def bt2_args = params.bowtie2_args ?: "--local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700"
    // Resolve index basename: if directory provided, find the index prefix
    """
    # Determine bowtie2 index prefix
    INDEX_PREFIX="${bowtie2_index_dir}"
    if [ -d "\${INDEX_PREFIX}" ]; then
        INDEX_PREFIX=\$(ls \${INDEX_PREFIX}/*.1.bt2 2>/dev/null | head -1 | sed 's/.1.bt2//')
        if [ -z "\${INDEX_PREFIX}" ]; then
            INDEX_PREFIX=\$(ls \${INDEX_PREFIX}/*.1.bt2l 2>/dev/null | head -1 | sed 's/.1.bt2l//')
        fi
    fi

    bowtie2 \\
        -p ${task.cpus} \\
        ${bt2_args} \\
        -x "\${INDEX_PREFIX}" \\
        -1 "${r1}" \\
        -2 "${r2}" \\
        2> "${sample_id}_bowtie2.log" \\
    | samtools view -bS -F 4 - \\
    > "${sample_id}.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*bowtie2-align-s version //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.bam ${sample_id}_bowtie2.log versions.yml
    """
}
