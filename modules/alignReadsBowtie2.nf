/*
 * Align reads using Bowtie2
 */

process alignReadsBowtie2 {
    
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }

    container 'ummidock/bowtie2_samtools:1.0.0-2'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path bt2Files  

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    echo "Running Bowtie2 Alignment for Reads"

    # Find index prefix from bt2 files
    INDEX=\$(ls *.1.bt2 | grep -v 'rev' | sed 's/\\.1\\.bt2\$//')
    echo "Index prefix: \$INDEX"

    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
        bowtie2 --no-unal -p ${task.cpus} -x "\$INDEX" \
        -1 "${reads[0]}" -2 "${reads[1]}" -S - \
        | samtools view -b - \
        | samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - \
        > "${sample_id}.bam"

    elif [ -f "${reads[0]}" ]; then
        bowtie2 --no-unal -p ${task.cpus} -x "\$INDEX" \
        -U "${reads[0]}" -S - \
        | samtools view -b - \
        | samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - \
        > "${sample_id}.bam"

    else
        echo "Error: Read file ${reads[0]} does not exist for sample ${sample_id}."
        exit 1
    fi

    echo "Alignment complete"
    """
}