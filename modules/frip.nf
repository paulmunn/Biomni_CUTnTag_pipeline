// ============================================================
// modules/frip.nf
// FRiP (Fraction of Reads in Peaks) calculation
// ============================================================
// FRiP = reads overlapping peaks / total mapped reads
// ENCODE recommends FRiP >= 0.01 for CUT&Tag; values > 0.2 are
// typical for high-quality H3K27ac or H3K4me3 experiments.

process FRIP_SCORE {
    tag "${sample_id}:${frip_type}"
    label 'process_low'

    publishDir "${params.outdir}/11_frip",
        mode: 'copy',
        pattern: "*.frip.txt"

    input:
    tuple val(sample_id), path(bam), path(bai), path(peak_file)
    val  frip_type

    output:
    tuple val(sample_id), path("${sample_id}_${frip_type}.frip.txt"), emit: tsv
    path "versions.yml",                                              emit: versions

    script:
    """
    # Total mapped reads (properly paired)
    TOTAL=\$(samtools view -c -F 4 -f 2 "${bam}")

    # Reads overlapping peaks
    IN_PEAKS=\$(bedtools intersect \\
        -a "${bam}" \\
        -b "${peak_file}" \\
        -u \\
        -f 0.20 \\
        | samtools view -c)

    # Calculate FRiP
    FRIP=\$(echo "scale=6; \${IN_PEAKS}/\${TOTAL}" | bc)

    echo -e "sample_id\\tfrip_type\\ttotal_reads\\treads_in_peaks\\tfrip_score" \\
        > "${sample_id}_${frip_type}.frip.txt"
    echo -e "${sample_id}\\t${frip_type}\\t\${TOTAL}\\t\${IN_PEAKS}\\t\${FRIP}" \\
        >> "${sample_id}_${frip_type}.frip.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "sample_id\\tfrip_type\\ttotal_reads\\treads_in_peaks\\tfrip_score" > ${sample_id}_${frip_type}.frip.txt
    echo -e "${sample_id}\\t${frip_type}\\t1000000\\t250000\\t0.25" >> ${sample_id}_${frip_type}.frip.txt
    touch versions.yml
    """
}
