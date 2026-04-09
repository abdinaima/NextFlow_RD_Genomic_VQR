process octopusCaller {
    
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }

    container 'dancooke/octopus:0.7.4'
    
    tag "$sample_id" 

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_fasta

    output:
    tuple val(sample_id), path("${sample_id}.vcf")  

    script:
    """
    echo "Running Octopus for sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "REF: ${genome_fasta}"

    octopus \\
        -R "${genome_fasta}" \\
        -I "${bam}" \\
        -o "${sample_id}.vcf"

    echo "Octopus complete: ${sample_id}.vcf"
    """
}

process indexOctopusVCF {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }

    container 'variantvalidator/gatk4:4.3.0.0'  

    tag "$sample_id"

    publishDir("$params.outdir/VCF", mode: "copy")  

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), path("${sample_id}.vcf.idx")  

    script:
    """
    echo "Indexing Octopus VCF for sample: ${sample_id}"

    gatk IndexFeatureFile -I "${sample_id}.vcf"

    echo "Indexing complete: ${sample_id}.vcf.idx"
    """
}