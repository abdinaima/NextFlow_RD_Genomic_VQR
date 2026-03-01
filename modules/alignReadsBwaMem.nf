/*
 * Align reads to the indexed genome
 */
process alignReadsBwaMem {   // Define the BWA-MEM process

    // Dynamically assigning a label to this process depending on which platform is specified in nextflow.config file. This will control which resource requirements are selected from the base.config file.
    if (params.platform == 'local') {                    // If workflow is run locally, then platform = "local" in nextflow.config file
        label 'process_low'                              // Process gets the label  'process_low' which means it will typically use few resources (e.g. CPU, memory, etc.)
    } else if (params.platform == 'cloud') {             // If workflow is run on the cloud, then platform = "cloud" in nextflow.config file
        label 'process_high'                             // Process gets the label  'process_high' which means more resources are assigned
    }
    container 'variantvalidator/indexgenome:1.1.0'   // Use the appropriate container with version explicitly stated

    tag "$sample_id"      // Tag jobs with sample ID for traceability

    input:
    tuple val(sample_id), path(reads)   // Input: 'reads' is a tuple of paths for paired-end reads
    path requiredIndexFiles             // All BWA MEM index files (.fasta.) needed for alignment. These files are staged into the process directory so BWA-MEM can use them to map reads to the reference genome.

    output:
    tuple val(sample_id), file("${sample_id}.bam")  // Output: this process will output a tuple: sample ID (value) and a single BAM named after the sample (file)

    script:
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')   # Finds the BWA-MEM index prefix by searching for a file ending in .amb and stripping the .amb extension. This gives the base name needed for my indexes, 

    echo "Running Align Reads"                                # Prints status message 
    echo "\$INDEX"                                            # Prints the index prefix being used

    # Check if the input FASTQ files exist
    if [ -f "${reads[0]}" ]; then                                         # Checks if first read file, reads[0], exists                 
        if [ -f "${reads[1]}" ]; then                                     # Then checks if also the second read file, reads[1], exists 
            # Paired-end                                                          # Runs BWA-MEM in paired-end mode
            bwa mem -M -t ${task.cpus} \$INDEX ${reads[0]} ${reads[1]} |          #BWA-MEM command: -M (marks shorter split hits as secondary), -t (number of threads/CPUs to use), INDEX (the base name of the BWA index), {reads[0]} {reads[1]} (the input FASTQ files for the forward and reverse reads) 
            samtools view -b - |                                                  # 'samtools view' is used to convert the SAM alignment file, to a BAM alignment file
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam         # 'samtools addreplacerg' command takes a BAM stream and adds read group metadara (e.g. read group ID, sample ID, Illlumina platform name), and then finally writes the results to the sample_id.bam. This can be important for downstream analysis tools that require read group information
        else                                                              
            # Single FASTQ mode                                            # But if only the first read file, reads[0], exists
            bwa mem -M -t ${task.cpus} \$INDEX ${reads[0]} |               # Then runs BWA-MEM in single-end mode
            samtools view -b - |                                           # 'samtools view' is used to convert the SAM alignment file, to a BAM alignment file
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam         # 'samtools addreplacerg' command takes a BAM stream and adds read group metadara (e.g. read group ID, sample ID, Illlumina platform name), and then finally writes the results to the sample_id.bam. This can be important for downstream analysis tools that require read group information
        fi
    else
        echo "Error: Read file ${reads[0]} does not exist for sample ${sample_id}."     # If reads[0] does not exist, meaning there are no reads at all, prints an error and exits.
        exit 1
    fi

    echo "Alignment complete"     # Prints a completion message.
    """
}
