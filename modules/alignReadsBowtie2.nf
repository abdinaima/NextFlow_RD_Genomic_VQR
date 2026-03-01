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
    tuple path(genomeFasta), path("index*.bt2")

    script:
    """
    echo "Running Bowtie2 Indexing"                       # Prints status message

    # Index the genome for Bowtie2 (output basename: index)
    bowtie2-build "${genomeFasta}" index

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
    path indexFiles                     // Input: Bowtie2 index files only
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")  // Output: this process will output a tuple: sample ID (value) and a single BAM named after the sample (file)

    
    script:
    """
    echo "Running Bowtie2 Alignment for Reads"                       # Prints status message 

    # Check if the input FASTQ files exist
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then                                # Checks if first read file, reads[0], exists  
        # Paired-end mode                                                    # Runs Bowtie2 in paired-end mode
        bowtie2 --no-unal -p ${task.cpus} -x index -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam \         #Bowtie2 command: --no-unal (optional argument, meaning reads that don't align to the reference genome will not be written to sam output), -p (number of threads/CPUs to use), -x (the genome index), -1 (input FASTQ file for the forward read), -2 (input FASTQ file for the reverse read), -S (output alignment in SAM format)
        samtools view -b - \                                                                                     # 'samtools view' is used to convert the SAM alignment file, to a BAM alignment file
        samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam     # 'samtools addreplacerg' command takes a BAM stream and adds read group metadara (e.g. read group ID, sample ID, Illlumina platform name), and then finally writes the results to the sample_id.bam. This can be important for downstream analysis tools that require read group information
    elif                                                                                                         # But if only the first read file, reads[0], exists
        # Single-end mode
        bowtie2 --no-unal -p ${task.cpus} -x index -U ${reads[0]} -S ${sample_id}.sam \     
        samtools view -b - \                                                                                     # 'samtools view' is used to convert the SAM alignment file, to a BAM alignment file
        samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam    # 'samtools addreplacerg' command takes a BAM stream and adds read group metadara (e.g. read group ID, sample ID, Illlumina platform name), and then finally writes the results to the sample_id.bam. This can be important for downstream analysis tools that require read group information
    
    else
        echo "Error: Read file ${reads[0]} does not exist for sample ${sample_id}."    # If reads[0] does not exist, meaning there are no reads at all, prints an error and exits.
        exit 1
    fi

    echo "Alignment complete"
    """
}