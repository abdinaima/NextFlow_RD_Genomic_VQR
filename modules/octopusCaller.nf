process octopusCaller {
    
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'dancooke/octopus:0.7.4'
    tag "$sample_id"
    shell = ['/bin/bash', '-c']  // ADD THIS LINE

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_all_files

    output:
    tuple val(sample_id), path("${sample_id}.vcf")

    script:
    """
    genome_fasta=\$(find . -maxdepth 1 -name '*.fasta' | head -1)
    echo "Running Octopus for sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "REF: \${genome_fasta}"
    ls -lh .

    octopus \
        -R "\${genome_fasta}" \
        --reads "${bam}" \
        -o "${sample_id}.vcf"

    echo "Octopus complete: ${sample_id}.vcf"
    """
}
