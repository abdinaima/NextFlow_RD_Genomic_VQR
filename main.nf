// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Print pipeline configuration
log.info """\
    ============================================
          DNASeq Pipeline Configuration
    ============================================
    platform        : ${params.platform}
    samplesheet     : ${params.samplesheet}
    genome          : ${params.genome_file}
    genome index    : ${params.genome_index_files}
    index genome    : ${params.index_genome}
    qsr truth vcfs  : ${params.qsrVcfs}
    output directory: ${params.outdir}
    fastp           : ${params.fastp}
    fastqc          : ${params.fastqc}
    aligner         : ${params.aligner}
    variant caller  : ${params.variant_caller}
    bqsr            : ${params.bqsr}
    degraded_dna    : ${params.degraded_dna}
    variant_recalibration: ${params.variant_recalibration}
    identity_analysis: ${params.identity_analysis}
    ============================================
""".stripIndent()

// Conditionally include modules
if (params.index_genome) {
    include { indexGenome } from './modules/indexGenome'
}
if (params.fastqc) {
    include { FASTQC } from './modules/FASTQC'
}
if (params.fastp) {
    include { fastp } from './modules/fastp'
}
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
if (params.bqsr) {
    include { baseRecalibrator } from './modules/BQSR'
}
include { combineGVCFs } from './modules/processGVCFs'
include { genotypeGVCFs } from './modules/processGVCFs'
if (params.variant_recalibration) {
    include { variantRecalibrator } from './modules/variantRecalibrator'
} else {
    include { filterVCF } from './modules/filterVCF'
}
if (params.identity_analysis) {
    include { identityAnalysis } from './modules/identityAnalysis'
}
if (params.aligner == 'bwa-mem') {
    include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
} else if (params.aligner == 'bwa-aln') {
    include { alignReadsBwaAln } from './modules/alignReadsBwaAln'
} else if (params.aligner == 'bowtie2') {
    include { alignReadsBowtie2 } from './modules/alignReadsBowtie2'
}
if (params.variant_caller == 'haplotype-caller') {
    include { haplotypeCaller } from './modules/haplotypeCaller'
} else {
    error "Unsupported variant caller: ${params.variant_caller}. Please specify 'haplotype-caller'."
}
if (params.degraded_dna) {
    include { mapDamage2 } from './modules/mapDamage'
    include { indexMapDamageBam } from './modules/indexBam'
}

workflow {

    // User decides to index genome or not
    if (params.index_genome){
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
    }
    else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
    }

    // Always create a separate fasta channel for processes that need it
    genome_fasta_ch = Channel.fromPath(params.genome_file)  // separate fasta channel

    // Create qsrc_vcf_ch channel
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)

    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            if (row.size() == 4) {
                tuple(row[0], [row[1], row[2]])
            } else if (row.size() == 3) {
                tuple(row[0], [row[1]])
            } else {
                error "Unexpected row format in samplesheet: $row"
            }
        }
    read_pairs_ch.view()

    // Run FASTQC on read pairs
    if (params.fastqc) {
        FASTQC(read_pairs_ch)
    }

    // Run fastp on read pairs and capture trimmed reads
    if (params.fastp) {
        fastp_out = fastp(read_pairs_ch)
        align_input_ch = fastp_out[0]  // trimmed reads passed forward
    } else {
        align_input_ch = read_pairs_ch  // fallback to untrimmed reads
    }

    // Align reads using trimmed or untrimmed reads depending on fastp param
    if (params.aligner == 'bwa-mem') {
        align_ch = alignReadsBwaMem(align_input_ch, indexed_genome_ch.collect())
    } else if (params.aligner == 'bwa-aln') {
        align_ch = alignReadsBwaAln(align_input_ch, indexed_genome_ch.collect())
    } else if (params.aligner == 'bowtie2') {
        bt2_index_ch = Channel.fromPath(params.genome_index_files).collect()
        align_ch = alignReadsBowtie2(align_input_ch, bt2_index_ch)
    }

    // Sort BAM files
    sort_ch = sortBam(align_ch)

    // Mark duplicates in BAM files
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files and collect the output channel
    indexed_bam_ch = indexBam(mark_ch)

    // Conditionally run mapDamage if degraded_dna parameter is set
    if (params.degraded_dna) {
        pre_mapDamage_ch = mapDamage2(indexed_bam_ch, indexed_genome_ch.collect())
        mapDamage_ch = indexMapDamageBam(pre_mapDamage_ch)
    } else {
        mapDamage_ch = indexed_bam_ch
    }

    // Create a channel from qsrVcfs
    knownSites_ch = Channel.fromPath(params.qsrVcfs)
        .filter { file -> file.getName().endsWith('.vcf.gz.tbi') || file.getName().endsWith('.vcf.idx') }
        .map { file -> "--known-sites " + file.getBaseName() }
        .collect()

    if (params.bqsr) {
        //  use genome_fasta_ch instead of indexed_genome_ch for BQSR
        bqsr_ch = baseRecalibrator(mapDamage_ch, knownSites_ch, genome_fasta_ch.collect(), qsrc_vcf_ch.collect())
    } else {
        bqsr_ch = mapDamage_ch
    }

    // Run HaplotypeCaller on BQSR files
    if (params.variant_caller == "haplotype-caller") {
        gvcf_ch = haplotypeCaller(bqsr_ch, genome_fasta_ch.collect()).collect()  // use genome_fasta_ch
    }

    // Now we map to create separate lists for sample IDs, VCF files, and index files
    all_gvcf_ch = gvcf_ch
        .collect { listOfTuples ->
            def sample_ids = listOfTuples.collate(3).collect { it[0] }
            def vcf_files = listOfTuples.collate(3).collect { it[1] }
            def vcf_index_files = listOfTuples.collate(3).collect { it[2] }
            return tuple(sample_ids, vcf_files, vcf_index_files)
        }

    // Combine GVCFs
    combined_gvcf_ch = combineGVCFs(all_gvcf_ch, genome_fasta_ch.collect())  // use genome_fasta_ch

    // Run GenotypeGVCFs
    final_vcf_ch = genotypeGVCFs(combined_gvcf_ch, genome_fasta_ch.collect())  // use genome_fasta_ch

    // Conditionally apply variant recalibration or filtering
    if (params.variant_recalibration) {
        def resourceOptions = [
            'Homo_sapiens_assembly38.known_indels': 'known=true,training=false,truth=false,prior=15.0',
            'hapmap_3.3.hg38': 'known=false,training=false,truth=true,prior=15.0',
            '1000G_omni2.5.hg38': 'known=false,training=true,truth=false,prior=12.0',
            '1000G_phase1.snps.high_confidence.hg38': 'known=true,training=true,truth=true,prior=10.0',
            'Homo_sapiens_assembly38.dbsnp138': 'known=true,training=false,truth=false,prior=2.0',
            'Mills_and_1000G_gold_standard.indels.hg38': 'known=true,training=true,truth=true,prior=12.0'
        ]
        knownSitesArgs_ch = Channel
            .fromPath(params.qsrVcfs)
            .filter { file -> file.getName().endsWith('.vcf.gz') || file.getName().endsWith('.vcf') }
            .map { file ->
                def baseName = file.getName().replaceAll(/\.vcf(\.gz)?$/, '')
                def resourceArgs = resourceOptions.get(baseName) ?: ""
                return "--resource:${baseName},${resourceArgs} ${file.getName()}"
            }
            .collect()
        filtered_vcf_ch = variantRecalibrator(final_vcf_ch, knownSitesArgs_ch, genome_fasta_ch.collect(), qsrc_vcf_ch.collect()) 
    } else {
        filtered_vcf_ch = filterVCF(final_vcf_ch, genome_fasta_ch.collect())  
    }

    // Conditionally run identityAnalysis if identity_analysis is true
    if (params.identity_analysis) {
        psam_info_ch = Channel
            .fromPath(params.samplesheet)
            .splitCsv(sep: '\t')
            .map { row ->
                if (row.size() == 4) {
                    tuple(row[0], row[3])
                } else if (row.size() == 3) {
                    tuple(row[0], row[2])
                } else {
                    error "Unexpected row format in samplesheet: $row"
                }
            }
        def combined_psam_content = new StringBuilder("#IID\tSID\tPAT\tMAT\tSEX\n")
        psam_file_ch = psam_info_ch.map { sample_info ->
            def sample_id = sample_info[0]
            def sex = sample_info[1]
            if (!sex) { sex = "NA" }
            def sample_line = "${sample_id}\t${sample_id}\t0\t0\t${sex}".stripIndent().trim()
            combined_psam_content.append(sample_line + "\n")
        }
        psam_file_ch.subscribe {
            def combined_psam_file = new File("/tmp/combined_samples.psam")
            combined_psam_file.text = combined_psam_content.toString()
            return combined_psam_file
        }
        identity_analysis_ch = identityAnalysis(filtered_vcf_ch, psam_file_ch)
    }
}

workflow FASTQC_only {
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            if (row.size() == 4) {
                tuple(row[0], [row[1], row[2]])
            } else if (row.size() == 3) {
                tuple(row[0], [row[1]])
            } else {
                error "Unexpected row format in samplesheet: $row"
            }
        }
    read_pairs_ch.view()
    if (params.fastqc) {
        FASTQC(read_pairs_ch)
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}