// ============================================================
// modules/pair_fastqs.nf
// Scan input directory, detect R1/R2 pairs, validate, output TSV
// ============================================================

process PAIR_FASTQS {
    tag "pair_fastqs"
    label 'process_single'

    publishDir "${params.outdir}/00_fastq_pairs", mode: 'copy'

    input:
    path input_dir
    val  paired_pattern

    output:
    path "fastq_pairs.tsv",   emit: pairs_tsv
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
