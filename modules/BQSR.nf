process baseRecalibrator {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'broadinstitute/gatk:4.1.4.0'
    tag "$bamFile"
    publishDir("$params.outdir/BAM", mode: "copy")
    
    input:
    tuple val(sample_id), file(bamFile), file(baiFile)
    val knownSites
    path indexFiles  // this receives genome_fasta_ch with fasta + fai + dict
    path qsrcVcfFiles

    output:
    tuple val(sample_id), file("${bamFile.baseName}_recalibrated.bam"), file("${bamFile.baseName}_recalibrated.bai")

    script:
    def knownSitesArgs = knownSites.join(' ')
    """
    echo "Running BQSR"

    # Find fasta using find instead of basename
    genomeFasta=\$(find -L . -name '*.fasta' | head -1)
    echo "Genome File: \${genomeFasta}"

    # Generate fai index if it doesn't exist
    if [ ! -f "\${genomeFasta}.fai" ]; then
        echo "Generating fasta index..."
        samtools faidx "\${genomeFasta}"
    fi

    # Generate dict if it doesn't exist
    if [ ! -f "\${genomeFasta%.*}.dict" ]; then
        echo "Generating sequence dictionary..."
        gatk CreateSequenceDictionary -R "\${genomeFasta}"
    fi

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    # Generate recalibration table for the input BAM file
    gatk --java-options "-Xmx8G" BaseRecalibrator \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        ${knownSitesArgs} \
        -O ${bamFile.baseName}.recal_data.table

    # Apply BQSR to the input BAM file
    gatk --java-options "-Xmx8G" ApplyBQSR \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        --bqsr-recal-file ${bamFile.baseName}.recal_data.table \
        -O ${bamFile.baseName}_recalibrated.bam

    # Index the recalibrated BAM file
    samtools index ${bamFile.baseName}_recalibrated.bam ${bamFile.baseName}_recalibrated.bai

    echo "BQSR Complete"
    """
}
