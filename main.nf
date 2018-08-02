params.reads = "$baseDir/data/samples/*.fastq"
params.genome = "$baseDir/data/genome.fa"
params.outdir = 'results'

log.info """\
         VARIANT CALLING TOY PIPELINE
         =============================
         genome: ${params.genome}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()

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
}
