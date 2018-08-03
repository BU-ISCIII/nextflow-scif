/*
========================================================================================
                  N E X T F L O W    B A S I C    W O R K F L O W
========================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/nextflow-scif
 @#### Authors
 S. Monzon <smonzon@isciii.es>

 ## Based on nf-core pipelines
 # https://github.com/nf-core
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:   BWA indexing of reference genome.
 - 2:   BWA mapping against indexed genome.
 - 3:   Samtools sorting and indexing.
 - 4:   Variant calling using bcftools
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     BU-ISCIII/nextflow-scif DEMO PIPELINE v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run BU-ISCIII/nextflow-scif -profile standard

    Pipeline arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --genome                  	Path to reference genome.
      --outdir						Output dir.
      --help						show this message.
      -profile                      Hardware config to use. standard/docker/singularity. Default: standard.

    """.stripIndent()
}

// Pipeline version
version = '0.1'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Default parameters
params.reads = "$baseDir/data/samples/*.fastq"
params.genome = "$baseDir/data/genome.fa"
params.outdir = 'results'


/*
 * the reference genome file
 */
genome_file = file(params.genome)

/*
 * Create the `reads` channel. Size 1 for single-end. Size 2 for paired-end.
 */
Channel
    .fromFilePairs( params.reads, size : 1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { reads }


// Header log info
log.info "========================================================="
log.info " BU-ISCIII/nextflow-scif : basic nf workflow"
log.info "========================================================="
def summary = [:]
summary['Reads']               = params.reads
summary['Reference genome']    = params.genome
summary['Results dir']         = params.outdir
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="


// Posible software version and profile checks (p.e check if standar profile is used in hpc server)
// if().....


/*
 * Step 1. Builds the genome index required by the mapping process
 */
process makeBWAindex {
	tag "${fasta}"

	input:
	file fasta from genome_file

	output:
	file "${fasta}*" into bwa_index

	script:
	"""
	mkdir BWAIndex
	bwa index -a bwtsw $fasta
	"""
}

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
    tag "$name"

    input:
    set val(name),file(reads) from reads
    file index from bwa_index
    file fasta from genome_file

    output:
    file '*.bam' into bwa_bam

    script:
    prefix = reads[0].toString() - ~/(\.fastq)$/
    """
    bwa mem -M $fasta $reads | samtools view -bS - > ${prefix}.bam
    """
}

process samtools {
    tag "${bam.baseName}"
    publishDir path: "${params.outdir}/bwa", mode: 'copy'

    input:
    file bam from bwa_bam

    output:
    file '*.sorted.bam' into bam_for_bcftools
    file '*.sorted.bam.bai' into bai_for_bcftools
    file '*.stats.txt' into samtools_stats

    script:
    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam -T ${bam.baseName}.sorted
    samtools index ${bam.baseName}.sorted.bam
    samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}.stats.txt
    """
}

process variantCalling {
	tag "${prefix}"
	publishDir path : "${params.outdir}/vcf", mode:'copy'

	input:
	file bam_sorted from bam_for_bcftools
	file bai_sorted from bai_for_bcftools
	file genome from genome_file
	output:
	file "*.vcf" into vcf_file

	script:
	prefix = bam_sorted[0].toString() - ~/(\.sorted)?(\.bam)?$/
	"""
	samtools mpileup -g -f ${genome} ${bam_sorted} | bcftools call -mv - > ${prefix}.vcf
	"""
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )

	// E-mail and html reporting configuration.

}
