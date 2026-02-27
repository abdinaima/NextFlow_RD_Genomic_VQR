/*
 * Run fastp trimming tool on the read fastq files
 */

process fastp { // Define the fastp process 

    container 'swglh/fastp:1.0.1' // Use the fastp container with version explicitly stated

    tag "$sample_id" // Tag jobs with sample ID for traceability

    publishDir("$params.outdir/fastp", mode: "copy")  // Output directory for fastp results

    input:
    tuple val(sample_id), path(reads) // Input: sample ID and reads (array of files). Reads is a tuple of paths for paired-end reads

    output:
    path "fastp_${sample_id}_logs/*"  // Output: This process will output all files to the path/directory that's named after the sample's log directory

    script:
    """
    echo "Running fastp"                                        # Print status message
    mkdir -p fastp_${sample_id}_logs                            # Create output directory for this sample

    # Check the number of files in reads and run fastp accordingly
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then                                # If both R1 and R2 exist (paired-end)
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R1.fastq.gz \
            -O fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R2.fastq.gz \
            -h fastp_${sample_id}_logs/fastp_${sample_id}.html \
            -j fastp_${sample_id}_logs/fastp_${sample_id}.json
    elif [ -f "${reads[0]}" ]; then                                                      # If only R1 exists (single-end)
        fastp -i ${reads[0]} \
            -o fastp_${sample_id}_logs/fastp_${sample_id}.trimmed.R1.fastq.gz \
            -h fastp_${sample_id}_logs/fastp_${sample_id}.html \
            -j fastp_${sample_id}_logs/fastp_${sample_id}.json
    else                                                                                 # If no valid reads found
        echo "No valid read files found for sample ${sample_id}"                         # Print error
        exit 1                                                                           # Exit with error
    fi

    echo "fastp Complete"                                                                # Print completion message
    """
}



