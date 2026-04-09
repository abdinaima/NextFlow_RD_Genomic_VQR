/*
 * Run fastp trimming tool on the read fastq files
 */

process fastp {

    container 'swglh/fastp:1.0.1'

    tag "$sample_id"

    publishDir("$params.outdir/fastp", mode: "copy")

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R1.fastq.gz"), path("fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R2.fastq.gz")  // trimmed reads passed forward
    path "fastp_${sample_id}_logs/*.{html,json}"  // logs published to outdir

    script:
    """
    echo "Running fastp"
    mkdir -p fastp_${sample_id}_logs

    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R1.fastq.gz \
            -O fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R2.fastq.gz \
            -h fastp_${sample_id}_logs/fastp_${sample_id}.html \
            -j fastp_${sample_id}_logs/fastp_${sample_id}.json
    elif [ -f "${reads[0]}" ]; then
        fastp -i ${reads[0]} \
            -o fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R1.fastq.gz \
            -h fastp_${sample_id}_logs/fastp_${sample_id}.html \
            -j fastp_${sample_id}_logs/fastp_${sample_id}.json
    else
        echo "No valid read files found for sample ${sample_id}"
        exit 1
    fi

    echo "fastp Complete"
    """
}