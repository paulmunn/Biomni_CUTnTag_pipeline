// ============================================================
// modules/pair_fastqs.nf
// Scan input directory (recursively) for paired-end FASTQ files,
// detect R1/R2 pairs, validate, and output a TSV manifest.
//
// IMPORTANT: input_dir is passed as `val` (a plain string), NOT
// as `path`. Using `path` would cause Nextflow to stage only the
// top-level directory contents into the work directory, breaking
// recursive subdirectory discovery. By passing the raw filesystem
// path as a string, pair_fastqs.py can call Path.rglob() freely
// across the full directory tree regardless of nesting depth.
// ============================================================

process PAIR_FASTQS {
    tag "pair_fastqs"
    label 'process_single'

    publishDir "${params.outdir}/00_fastq_pairs", mode: 'copy'

    input:
    val  input_dir       // Raw filesystem path string — NOT staged as path
    val  paired_pattern

    output:
    path "fastq_pairs.tsv",    emit: pairs_tsv
    path "pairing_report.txt", emit: report

    script:
    """
    pair_fastqs.py \\
        --input_dir "${input_dir}" \\
        --pattern   "${paired_pattern}" \\
        --output    fastq_pairs.tsv \\
        --report    pairing_report.txt
    """

    stub:
    """
    touch fastq_pairs.tsv pairing_report.txt
    """
}
