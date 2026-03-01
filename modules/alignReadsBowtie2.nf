/*
 * Align reads using Bowtie2
 */

process indexGenomeBowtie2 {     // Defines the Bowtie2 indexing process
    
    // Dynamically assigning a label to this process depending on which platform is specified in nextflow.config file. This will control which resource requirements are selected from the base.config file.
    if (params.platform == 'local') {              // If workflow is run locally, then platform = "local" in nextflow.config file
        label 'process_low'                        // Process gets the label  'process_low' which means it will typically use few resources (e.g. CPU, memory, etc.)
    } else if (params.platform == 'cloud') {       // If workflow is run on the cloud, then platform = "cloud" in nextflow.config file
        label 'process_high'                       // Process gets the label  'process_high' which means more resources are applied
    }

    container 'ummidock/bowtie2_samtools:1.0.0-2' // Use the appropriate container with version explicitly stated

    tag "${genomeFasta.baseName}" // Tag jobs with sample ID for traceability

    input:
    path genomeFasta

    output:
    tuple val("bt2_index"), path("bt2_index*.bt2")

    script:
    """
    echo "Running Bowtie2 Indexing"                       # Prints status message

    # Index the genome for Bowtie2
    bowtie2-build "${genomeFasta}" bt2_index
    ls -lh bt2_index*.bt2

    echo "Bowtie2 indexing complete"

    """
}



process alignReadsBowtie2 {   // Defines the Bowtie2 alignment process
   
    // Dynamically assigning a label to this process depending on which platform is specified in nextflow.config file. This will control which resource requirements are selected from the base.config file.
    if (params.platform == 'local') {              // If workflow is run locally, then platform = "local" in nextflow.config file
        label 'process_low'                        // Process gets the label  'process_low' which means it will typically use few resources (e.g. CPU, memory, etc.)
    } else if (params.platform == 'cloud') {       // If workflow is run on the cloud, then platform = "cloud" in nextflow.config file
        label 'process_high'                       // Process gets the label  'process_high' which means more resources are applied
    }

    container 'ummidock/bowtie2_samtools:1.0.0-2' // Use the appropriate container with version explicitly stated

    tag "$sample_id" // Tag jobs with sample ID for traceability

    input:
    tuple val(sample_id), path(reads)   // Input: 'reads' is a tuple of paths for paired-end reads
    tuple val(indexPrefix), path(bt2Files)                     // Input: Bowtie2 index files only
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")  // Output: this process will output a tuple: sample ID (value) and a single BAM named after the sample (file)

    
    script:
    """
    echo "Running Bowtie2 Alignment for Reads"                       # Prints status message 
    echo "Index prefix: ${indexPrefix}"

    # Check if the input FASTQ files exist
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
        # paired-end
        bowtie2 --no-unal -p ${task.cpus} -x "${indexPrefix}" \
        -1 "${reads[0]}" -2 "${reads[1]}" -S - \
        | samtools view -b - \
        | samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - \
        > "${sample_id}.bam"

    elif [ -f "${reads[0]}" ]; then
        # single-end
        bowtie2 --no-unal -p ${task.cpus} -x "${indexPrefix}" \
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