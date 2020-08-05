#!/usr/bin/env nextflow

/*****************************************************  
 *--------------------NEXTFLOW-----------------------*
 *****************************************************  
 * NGS pipeline processes for GATK CNV analysis      * 
 * Can be configured using nextflow settings         * 
 ***************************************************** 
 * Version : v1.0.0                                  * 
 ***************************************************** 
 * #### Authors                                      *
 * Benjamin NAVET     							     * 
 *****************************************************/

HOME="/media/WDISK1"
srcFolder="/CDISK2/TEST-CNV-BOUR"

EMPLACEMENT_TOOLS="$HOME/tools"
EMPLACEMENT_DB="$HOME/db"
EMPLACEMENT_SCRIPTS="$HOME/scripts"


PICARD_TOOL="$EMPLACEMENT_TOOLS/picard.jar"
GATK="$EMPLACEMENT_TOOLS/gatk-4.1.5.0/gatk"
BWA="$EMPLACEMENT_TOOLS/bwa-0.7.17/bwa"

REFERENCE="$EMPLACEMENT_DB/HG38_outalt/Homo_sapiens_assembly38_only_autochr.fasta"
KNOWN_SITES="$EMPLACEMENT_DB/HG38_GATK/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
TARGET="$EMPLACEMENT_DB/S30409818_Covered.bed"

params.bams	= "/CDISK2/TEST-CNV-BOUR/*.{bam,bam.bai}"
bam_ch		= Channel.fromFilePairs(params.bams)

process PreprocessIntervals {

	output:
		file("targets.preprocessed.interval_list") into (Kit_ch)
	script: 
	"""
	$GATK PreprocessIntervals -R $REFERENCE -L $TARGET --bin-length 0 -imr OVERLAPPING_ONLY -O targets.preprocessed.interval_list
	"""


}

process AnnotateIntervals {
	publishDir "${srcFolder}/READCOUNT/", mode : 'copy'	
	input:
		file("targets.preprocessed.interval_list") from Kit_ch
	output: 
		file("targets.annotated.tsv") into (annoKit_ch)
	script:
	"""
	$GATK AnnotateIntervals -L targets.preprocessed.interval_list -R $REFERENCE -imr OVERLAPPING_ONLY -O targets.annotated.tsv
	"""
}


process CollectReadCounts {
	tag "${sampleId}"
    	publishDir "${srcFolder}/READCOUNT/", mode : 'copy'
	maxForks 32
	cpus 1

	input:
		set sampleId, file(bam) from bam_ch
		file("targets.preprocessed.interval_list") from Kit_ch
	output:	file "${sampleId}.hdf5" into read_ch, read_ch_Deter
	script:
	"""
	$GATK CollectReadCounts -L targets.preprocessed.interval_list -R $REFERENCE -imr OVERLAPPING_ONLY -I ${bam[0]} -O ${sampleId}.hdf5 
	"""

}


process FilterIntervals {
	publishDir "${srcFolder}/READCOUNT/", mode : 'copy'	
	output:
		file("cohort.gc.filtered.interval_list") into (cohorte_gc)
	input: 
		file("targets.preprocessed.interval_list") from Kit_ch
		file("targets.annotated.tsv") from annoKit_ch
	script:
	"""
	$GATK FilterIntervals -L targets.preprocessed.interval_list --annotated-intervals targets.annotated.tsv -imr OVERLAPPING_ONLY --output cohort.gc.filtered.interval_list
	"""

}

process DetermineGermlineContigPloidy {
	conda '/WDISK1/anaconda3/envs/gatk'
	publishDir "${srcFolder}/MODEL/", mode : 'copy'	
	cpus 10
	input: 
		file(read) from read_ch_Deter.collect()
		file("cohort.gc.filtered.interval_list") from cohorte_gc

	output:
		file 'ploidy*' into ch_ploidy
	script:
	"""
	$GATK DetermineGermlineContigPloidy -L cohort.gc.filtered.interval_list --interval-merging-rule OVERLAPPING_ONLY -I ${read.join(' -I ')} --contig-ploidy-priors $EMPLACEMENT_DB/table_ploidy_autosex.tsv --output . --output-prefix ploidy
	"""

}

process GermlineCNVCaller {
	conda '/WDISK1/anaconda3/envs/gatk'
	publishDir "${srcFolder}/RESULT/", mode : 'copy'	
	input: 
		file("targets.preprocessed.interval_list") from Kit_ch
		file(read) from read_ch.collect()
		file (ploidy) from ch_ploidy
		file("targets.annotated.tsv") from annoKit_ch

	output:
	file 'cohort24*' into ch_ploidy2

	script:
	"""
	$GATK GermlineCNVCaller --run-mode COHORT -L targets.preprocessed.interval_list --interval-merging-rule OVERLAPPING_ONLY \
	-I ${read.join(' -I ')} --contig-ploidy-calls ${srcFolder}/MODEL/ploidy-calls --annotated-intervals targets.annotated.tsv --output . --output-prefix cohort24 --verbosity DEBUG
	"""

}


/*Faire un sous nextflow
process PostprocessGermlineCNVCalls {
	conda '/WDISK1/anaconda3/envs/gatk'
	publishDir "${srcFolder}/RESULT_SAMPLE/", mode : 'copy'	

	output:
	file '*' into ch_ploidy25
	
	script:
	"""
	$GATK PostprocessGermlineCNVCalls --model-shard-path ${srcFolder}/RESULT/cohort24-model --calls-shard-path ${srcFolder}/RESULT/cohort24-calls --allosomal-contig chrX --allosomal-contig chrY \
--contig-ploidy-calls ${srcFolder}/MODEL/ploidy-calls --sample-index 2 --output-genotyped-intervals Sample_interval_TEST.vcf.gz --output-genotyped-segments Sample_segments_TEST.vcf.gz --output-denoised-copy-ratios Sample_denoidsed_TEST --sequence-dictionary $EMPLACEMENT_DB/HG38_outalt/Homo_sapiens_assembly38_only_autochr.dict
	"""


}

*/

