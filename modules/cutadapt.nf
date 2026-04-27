// ============================================================
// modules/cutadapt.nf
// Trim adapters from paired-end reads using Cutadapt
// ============================================================
// Recommendation: Nextera transposase adapters are the default
// for CUT&Tag (CTGTCTCTTATACACATCT). Override with --adapter_fwd
// and --adapter_rev if using a different library prep.

process CUTADAPT {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/02_trimmed",
        mode: 'copy',
        pattern: "*.fastq.gz"
    publishDir "${params.outdir}/02_trimmed/logs",
        mode: 'copy',
        pattern: "*.log"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_cutadapt.log",  emit: log
    path "versions.yml",               emit: versions

    script:
    def adapter_fwd = params.adapter_fwd ?: "CTGTCTCTTATACACATCT"
    def adapter_rev = params.adapter_rev ?: "CTGTCTCTTATACACATCT"
    def min_len     = params.min_length  ?: 20
    """
    cutadapt \\
        --cores ${task.cpus} \\
        -a "${adapter_fwd}" \\
        -A "${adapter_rev}" \\
        --minimum-length ${min_len} \\
        --quality-cutoff 20 \\
        --trim-n \\
        --discard-trimmed \\
        -o "${sample_id}_R1_trimmed.fastq.gz" \\
        -p "${sample_id}_R2_trimmed.fastq.gz" \\
        "${r1}" "${r2}" \\
        > "${sample_id}_cutadapt.log" 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_R1_trimmed.fastq.gz ${sample_id}_R2_trimmed.fastq.gz
    touch ${sample_id}_cutadapt.log versions.yml
    """
}
