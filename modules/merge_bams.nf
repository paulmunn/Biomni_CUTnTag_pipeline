// ============================================================
// modules/merge_bams.nf
// Merge BAM files for group-level peak calling
// ============================================================

process MERGE_BAMS {
    tag "${group_id}"
    label 'process_medium'

    publishDir "${params.outdir}/09_peaks_merged_groups/${group_id}",
        mode: 'copy',
        pattern: "*.{bam,bai}"

    input:
    tuple val(group_id), path(bams), path(bais), val(bam_type)

    output:
    tuple val(group_id), path("${group_id}_merged.bam"), path("${group_id}_merged.bam.bai"), val(bam_type), emit: merged_bam
    path "versions.yml", emit: versions

    script:
    def bam_list = bams instanceof List ? bams.join(' ') : bams
    def n_bams   = bams instanceof List ? bams.size() : 1
    """
    if [ "${n_bams}" -eq 1 ]; then
        # Only one BAM: just copy and index
        cp "${bam_list}" "${group_id}_merged.bam"
    else
        samtools merge \\
            -@ ${task.cpus} \\
            -f \\
            "${group_id}_merged.bam" \\
            ${bam_list}
    fi

    samtools sort \\
        -@ ${task.cpus} \\
        -m 2G \\
        -o "${group_id}_merged_sorted.bam" \\
        "${group_id}_merged.bam"

    mv "${group_id}_merged_sorted.bam" "${group_id}_merged.bam"

    samtools index \\
        -@ ${task.cpus} \\
        "${group_id}_merged.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${group_id}_merged.bam ${group_id}_merged.bam.bai versions.yml
    """
}
