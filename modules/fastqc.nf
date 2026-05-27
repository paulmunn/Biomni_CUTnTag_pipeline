// ============================================================
// modules/fastqc.nf
// Run FastQC on paired-end FASTQ files
// ============================================================
// Output design note:
//   FastQC outputs are declared with EXPLICIT file names (no globs).
//   Glob outputs (path("*.zip")) emit DataflowStream objects in DSL2,
//   which carry a PoisonPill terminator that propagates into .collect()
//   and causes "Not a valid path value type: DataflowStreamReadAdapter".
//   By renaming inputs to known names and declaring exact output paths,
//   we guarantee plain Path values are emitted.
// ============================================================

process FASTQC {
    tag "${sample_id}:${stage}"
    label 'process_low'

    publishDir "${params.outdir}/0${stage == 'raw' ? '1' : '3'}_fastqc_${stage}",
        mode: 'copy',
        pattern: "*.{html,zip}"

    input:
    tuple val(sample_id), path(r1), path(r2)
    val  stage

    output:
    // Explicit named outputs — NO globs — to avoid DataflowStream emission
    tuple val(sample_id), path("${sample_id}_${stage}_R1_fastqc.html"),
                          path("${sample_id}_${stage}_R2_fastqc.html"), emit: html
    tuple val(sample_id), path("${sample_id}_${stage}_R1_fastqc.zip"),
                          path("${sample_id}_${stage}_R2_fastqc.zip"),  emit: zip
    path "versions.yml",                                                 emit: versions

    script:
    def prefix = "${sample_id}_${stage}"
    """
    # Rename inputs to include stage tag for MultiQC disambiguation
    ln -sf "${r1}" "${prefix}_R1.fastq.gz"
    ln -sf "${r2}" "${prefix}_R2.fastq.gz"

    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        "${prefix}_R1.fastq.gz" \\
        "${prefix}_R2.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_${stage}_R1_fastqc.html ${sample_id}_${stage}_R1_fastqc.zip
    touch ${sample_id}_${stage}_R2_fastqc.html ${sample_id}_${stage}_R2_fastqc.zip
    touch versions.yml
    """
}
