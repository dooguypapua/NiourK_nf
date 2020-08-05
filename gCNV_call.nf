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

resultFolter="/CDISK2"

cohorte_model_determine_ploidy="/media/CDISK2/TEST-CNV-BOUR/MODEL/ploidy-model/"
cohorte_model_gcnv_caller="/media/CDISK2/TEST-CNV-BOUR/RESULT/cohort24-model"

/*Chemin*/
EMPLACEMENT_TOOLS="$HOME/tools"
EMPLACEMENT_DB="$HOME/db"
EMPLACEMENT_SCRIPTS="$HOME/scripts"

/*Outils*/
GATK="$EMPLACEMENT_TOOLS/gatk-4.1.5.0/gatk"

REFERENCE="$EMPLACEMENT_DB/HG38_outalt/Homo_sapiens_assembly38_only_autochr.fasta"
TARGET="$EMPLACEMENT_DB/S30409818_Covered.bed"

params.bams	= "/media/CDISK2/ANALYSE_ESTELLE/check-sex-master/*.{bam,bam.bai}"
bam_ch		= Channel.fromFilePairs(params.bams)


/* creation du ficier d'intervalle a partir du fichier bed*/

process PreprocessIntervals {
	output:
		file("targets.preprocessed.interval_list") into (Kit_ch)
	script: 
	"""
	$GATK PreprocessIntervals -R $REFERENCE -L $TARGET --bin-length 0 --padding 250 -imr OVERLAPPING_ONLY -O targets.preprocessed.interval_list
	"""
}


/*
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
*/

/* Compte les reads pour chaque echantillon */
process CollectReadCounts {
	tag "${sampleId}"
    	publishDir "${resultFolter}/READCOUNT/", mode : 'copy'
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


/*CASE MODE */

process DetermineGermlineContigPloidy {
	conda '/WDISK1/anaconda3/envs/gatk'
	publishDir "${resultFolter}/RESULT_CASE_MODE/CALL", mode : 'copy'	
	cpus 10
	input: 
		file(read) from read_ch_Deter

	output:
		file 'ploidy-case*' into ch_ploidy
	script:
	"""
	$GATK DetermineGermlineContigPloidy --model ${cohorte_model_determine_ploidy} -I ${read} --output ${read.baseName} --output-prefix ploidy-case --verbosity DEBUG
	"""

}

/*
process GermlineCNVCaller {
	conda '/WDISK1/anaconda3/envs/gatk'
	publishDir "${resultFolter}/RESULT_CASE_MODE/RESULT", mode : 'copy'	
	input: 
		file(read) from read_ch
		file (ploidy) from ch_ploidy

	output:
	file '${read.baseName}' into ch_ploidy2

	script:
	"""
	$GATK GermlineCNVCaller --run-mode CASE -I ${read} --contig-ploidy-calls ${resultFolter}/RESULT_CASE_MODE/CALL/ploidy-case-calls/ --model ${cohorte_model_gcnv_caller} --output ${read.baseName} --output-prefix case --verbosity DEBUG
	"""

}
*/
/*
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

