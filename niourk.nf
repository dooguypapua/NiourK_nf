#!/usr/bin/env nextflow


/*
**********************************************************************
                             HELP & USAGE                             
**********************************************************************
*/
// Colors
c_reset  = params.monochrome_logs ? '' : "\033[0m";
c_dim    = params.monochrome_logs ? '' : "\033[2m";
c_title = params.monochrome_logs ? '' : "\033[38;2;71;157;230m";
c_subtitle = params.monochrome_logs ? '' : "\033[38;2;75;138;227m";
c_light1 = params.monochrome_logs ? '' : "\033[38;2;252;146;86m";
c_light2 = params.monochrome_logs ? '' : "\033[38;2;200;200;200m";
c_grey = params.monochrome_logs ? '' : "\033[38;2;150;150;150m";
c_red = params.monochrome_logs ? '' : "\033[38;2;212;64;89m";
c_green = params.monochrome_logs ? '' : "\033[38;2;55;200;113m";

def helpMessage(strError) {
    log.info niourkHeader()
    println c_reset
    log.error"""
    $strError
    """
    log.info"""
    .
    [Usage]

      The typical command for running the pipeline is as follows:

      nextflow run niourk.nf --pathBam bamfolder --pathOut results --genome hg38
      
      Mandatory arguments:
        --path_bam                  Path to input BAMs folder
        --path_out                  Path to results folder
        --genome                    Reference Genome (hg19, hg38, rCRS, ...)
        --mito                      Mitochondrial mode (true|false)

      Sequencing arguments:
        --path_bed                  Path to target BED file
        --path_param                Path to sequencing parameters files
                                    (Torrent  = ion_params_00.json)
                                    (Illumina = RunParameters.xml )

      Tools arguments:
        --path_gatk                 Path to Genome Analysis Toolkit (GATK) executable
        --path_tvc                  Path to Torrent Variant Caller (TVC) executable
        --path_strelka              Path to Strelka executable
        --path_samtools             Path to samtools executable
        --path_bedtools             Path to bedtools executable
        --path_vt                   Path to vt executable
        --path_vep                  Path to VEP executable
        --path_vep_cache            Path to VEP cache directory
        --path_vcfanno              Path to VcfAnno executable
        --path_vcfvalidator         Path to EBIvariation vcf-validator
        --version_deepvariant       Version tag for deepvariant docker container

      Calling arguments:
        --min_baseq                 Minimum base quality to consider a call
        --min_mapq                  Minimum mapping quality to consider a call
        --min_cov                   Minimum coverage to consider a call
        --min_af                    Minimum variant allele frequency
        --max_sb                    Maximum variant strand-bias
        --min_varcov                Mininimum variant coverage
        --min_varscore              Minimum variant score

      Other options:
        --help                      Print help
        --path_log                  Path to log folder
        --monochrome_logs           Logs will be without colors
    .
    """.stripIndent()
}
// Show help message
if (params.help) exit 0, helpMessage("")

// Check input & output
if (!params.path_bam) exit 1, helpMessage(".\n[ERROR]\n  Parameter `--path_bam` is required")
bamFiles = file(params.path_bam).list()
bamCount = 0
for( def bam : bamFiles ) { if (file(bam).getExtension()=="bam") bamCount+=1 }
if (!params.path_out) exit 1, helpMessage(".\n[ERROR]\n  Parameter `--path_out` is required")
pathVCF = file(params.path_out+"/VCFs")
if (!pathVCF.isDirectory()) pathVCF.mkdir()
pathVCFraw = file(params.path_out+"/VCFs/raw")
if (!pathVCFraw.isDirectory()) pathVCFraw.mkdir()

// Check reference
if (params.genome=="") exit 1, helpMessage(".\n[ERROR]\n  Please specify a reference `--genome`")
pathGenomeDir = "${workflow.projectDir}/reference/${params.genome}"
if (!file(pathGenomeDir).isDirectory()) exit 1, helpMessage(".\n[ERROR]\n  Invalid genome `"+workflow.projectDir+"/reference/"+params.genome+"` not found")
pathFasta = file(pathGenomeDir+"/"+params.genome+".fasta")
pathFastaFai = file(pathGenomeDir+"/"+params.genome+".fasta.fai")
if (!pathFasta.exists()) exit 1, helpMessage(".\n[ERROR]\n  Reference fasta not found `"+params.path_fasta+"` not found")
if (!pathFastaFai.exists()) exit 1, helpMessage(".\n[ERROR]\n  Reference fasta not found `"+params.path_fastaFai+"` not found")

// Mitochondrial mode
mitoMode = params.mito.toBoolean()
if ((params.mito!=false) && (params.mito!=true)) exit 1, helpMessage(".\n[ERROR]\n  Invalid mitochondrial mode `"+params.mito+"`, specify `true` or `false`")

// Check sequencing parameters (if not specified, search in bam folder)
pathParam = params.path_param
platform = null
if (pathParam==null) {
  if (file("${params.path_bam}/ion_params_00.json").exists()) pathParam = "${params.path_bam}/ion_params_00.json" ; 
  if (file("${params.path_bam}/RunParameters.xml").exists()) pathParam = "${params.path_bam}/RunParameters.xml"
}
if ((pathParam!=null) && (!file(pathParam).exists())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid path for `--path_param`")
nameParam = file(pathParam).getName()
// Copy parameters file to outdir
if (pathParam!=null) file(pathParam).copyTo("${params.path_out}/${nameParam}")
// Platform affectation
if (pathParam.contains("json")) platform = "IonTorrent"
else platform = "Illumina"
// Target BED file
//@ must be sorted + one version bgzip -f ${pathBed} && tabix -p bed ${pathBed}.gz
pathBed = params.path_bed
targetName = file(pathBed).getName()
if ((pathBed!=null) && (!file(pathBed).exists())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid path for `--path_bed`")
if (file(pathBed).getExtension()=="bed") pathBedgz = formatTargetBed()

// Check tools path
missing_tools = ""
if ((params.path_gatk=="") || (!file(params.path_gatk).exists())) missing_tools+="\n    `--path_gatk`"
if ((params.path_tvc=="") || (!file(params.path_tvc).exists())) missing_tools+="\n   `--path_tvc`"
if ((params.path_strelka=="") || (!file(params.path_strelka).exists())) missing_tools+="\n   `--path_strelka`"
if ((params.path_samtools=="") || (!file(params.path_samtools).exists())) missing_tools+="\n   `--path_samtools`"
if ((params.path_bedtools=="") || (!file(params.path_bedtools).exists())) missing_tools+="\n   `--path_bedtools`"
if ((params.path_vt=="") || (!file(params.path_vt).exists())) missing_tools+="\n   `--path_vt`"
if ((params.path_vep=="") || (!file(params.path_vep).exists())) missing_tools+="\n   `--path_vep`"
if ((params.path_vep_cache=="") || (!file(params.path_vep_cache).isDirectory())) missing_tools+="\n   `--path_vep_cache`"
if ((params.path_vcfanno=="") || (!file(params.path_vcfanno).exists())) missing_tools+="\n   `--path_vcfanno`"
if ((params.path_vcfvalidator=="") || (!file(params.path_vcfvalidator).exists())) missing_tools+="\n   `--path_vcfvalidator`"
if (missing_tools!="") exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid path for params: "+missing_tools)
if ((params.version_deepvariant=="") || (!params.version_deepvariant)) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid version for params: `--version_deepvariant`")

// Calling
if ((params.min_baseq=="") || (!params.min_baseq.toString().isNumber())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid integer for `--min_baseq`")
if ((params.min_mapq=="") || (!params.min_mapq.toString().isNumber())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid integer for `--min_mapq`")
if ((params.min_af=="") || (!params.min_af.toString().isFloat())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid float for `--min_af`")
if ((params.min_cov=="") || (!params.min_cov.toString().isNumber())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid integer for `--min_cov`")
if ((params.min_varcov=="") || (!params.min_varcov.toString().isNumber())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid integer for `--min_varcov`")
if ((params.min_varscore=="") || (!params.min_varscore.toString().isNumber())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid integer for `--min_varscore`")
if ((params.max_sb=="") || (!params.max_sb.toString().isFloat())) exit 1, helpMessage(".\n[ERROR]\n  Please specify a valid float for `--max_sb`")


// VEP Cache directory
if (!file(params.path_vep_cache+"/homo_sapiens_refseq").isDirectory())  exit 1, helpMessage(".\n[ERROR]\n  Any VEP cache directory found `"+params.path_vep_cache+"/homo_sapiens_refseq`")
if (!file(params.path_vep_cache+"/homo_sapiens_refseq/100_"+params.genome).isDirectory())  exit 1, helpMessage(".\n[ERROR]\n  Any VEP cache directory found `"+params.path_vep_cache+"/homo_sapiens_refseq/100_"+params.genome+"`")

// Channels
Channel
    .fromFilePairs(params.path_bam+'/*.{bam,bai}') { file -> file.name.replaceAll(/.bam|.bai$/,'') }
    .set { ch_Sample }

// split bam channel into multiples
ch_Sample.into { ch_SampleToValidate ; ch_SampleToGATKhc ; ch_SampleToMutect2 ; ch_SampleToTVC ; ch_SampleToDeepvariant ; ch_SampleToStrelka2 ; ch_SampleVCF }


test = true
/*
**********************************************************************
                           DISPLAY & SUMMARY                          
**********************************************************************
*/
// Summary
def summary = [:]
summary['User']       = workflow.userName
summary['nfDir'] = workflow.launchDir
summary['Input']      = params.path_bam
summary['Output']     = params.path_out
summary['Genome']     = params.genome+"\n― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ― ―"
summary['BAMs']       = bamCount+" file(s) found"
summary['Platform']   = platform
summary['Param']      = pathParam
summary['Target']     = pathBed
log.info niourkHeader()
print c_reset
if (!params.monochrome_logs) print c_light1
println "―――――――――――――――――――――――――――――――――――――――――――――――――――――――――"
log.info summary.collect { k, v -> "${k.padRight(10)}: $v" }.join("\n")
println "―――――――――――――――――――――――――――――――――――――――――――――――――――――――――${c_reset}"

// Niourk header
def niourkHeader() {""
    //printf "${params.color_title} ███╗   ██╗██╗ ██████╗ ██╗  \e[0m"
    println c_title
    return """
███╗   ██╗██╗ ██████╗ ██╗   ██╗██████╗ ██╗  ██╗
████╗  ██║██║██╔═══██╗██║   ██║██╔══██╗██║ ██╔╝
██╔██╗ ██║██║██║   ██║██║   ██║██████╔╝█████╔╝ 
██║╚██╗██║██║██║   ██║██║   ██║██╔══██╗██╔═██╗ 
██║ ╚████║██║╚██████╔╝╚██████╔╝██║  ██║██║  ██╗
╚═╝  ╚═══╝╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝ v.${workflow.manifest.version}
    """.stripIndent()
}



/*
**********************************************************************
                       CHECKING REQUIRED FILES                        
**********************************************************************
*/
def formatTargetBed() {
    cmdBed = "bedtools sort -i "+pathBed+" > work/"+file(pathBed).getBaseName()+"_sort.bed && bgzip -f work/"+file(pathBed).getBaseName()+"_sort.bed && tabix -p bed work/"+file(pathBed).getBaseName()+"_sort.bed.gz"
    def process = [ 'bash', '-c', cmdBed ].execute()
    process.waitFor()
    if (process.exitValue()!=0) exit 1, helpMessage(".\n[ERROR]\n  Unable to sort/bgzip/index `"+pathBed+"`")
    return "work/"+file(pathBed).getBaseName()+"_sort.bed.gz"
}



/* Parse software version numbers */
process GetSoftwareVersions {
    // publishDir The publishDir directive allows you to publish the process output files to a specified folder
    script:
      """
      python3 $workflow.projectDir/scripts/Nk_checkTools.py gatk=${params.path_gatk} tvc=${params.path_tvc} \
        strelka=${params.path_strelka} samtools=${params.path_samtools} bedtools=${params.path_bedtools} \
        vt=${params.path_vt} vcfanno=${params.path_vcfanno} vcf-validator=${params.path_vcfvalidator} \
        deepvariant=${params.version_deepvariant} > ${params.path_out}/tools_version.json
      """
}


/* Validate & index input BAM files */
process ValidateBamFile {
    label 'cpus_8' // allow group process with same computing reuqirements

    tag "$sample" // tag in report tasks

    errorStrategy "terminate"

    input:
      set sample, file(bam) from ch_SampleToValidate

    output:
      file("${sample}.bamidxstats") into ch_BamIndexStatsPlot

    script:
      java_options = "\"-Xms"+(task.memory.toGiga() / 2).trunc()+"g -Xmx"+(task.memory.toGiga())+"g\""
      """
      ${params.path_gatk} --java-options ${java_options} ValidateSamFile -I ${sample}.bam >/dev/null 2>&1 || (>&2 echo "[ERROR] ValidateSamFile" && exit 1)


      python3 $workflow.projectDir/scripts/Nk_checkBam.py refCompatibility ${params.path_samtools} ${task.cpus} ${sample}.bam $workflow.projectDir/reference/${params.genome}/${params.genome}.dict

      ${params.path_gatk} --java-options ${java_options} BamIndexStats -I ${sample}.bam > ${sample}.bamidxstats

      """
}



/*      
**********************************************************************
                           BAM MANIPULATION                           
**********************************************************************
*/
/* Index and BamIndexStats plot */
process BamIndexStatsPlot {
    label 'cpus_8' // allow group process with same computing reuqirements

    tag "$bamidxstats.baseName" // tag in report tasks

    input:
      file bamidxstats from ch_BamIndexStatsPlot.collect()

    output:
      file ("BamValidated.txt") into ch_BamValidated

    script:
      """
      tail -n +1 ${bamidxstats} > ${params.path_out}/bamIdxStats.txt
      python3 $workflow.projectDir/scripts/Nk_checkBam.py plotBamidxstats ${params.path_out}/bamIdxStats.txt
      echo "SUCCESS" > BamValidated.txt
      """
}



/*      
**********************************************************************
                               CALLING                                
**********************************************************************
// All VCFs are created in their work dir before validation and copy to pathVCFraw
*/
/* GATK HaplotypeCaller */
process HaplotypeCaller {
    label 'cpus_max' // allow group process with same computing reuqirements

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam) from ch_SampleToGATKhc
      file bamvalidated from ch_BamValidated

    output:
      file("${sample}_gatkHC.vcf") into ch_vcfGATKhc

    when: (test==true) //mitoMode==false

    script:
    java_options = "\"-Xms"+(task.memory.toGiga() / 2).trunc()+"g -Xmx"+(task.memory.toGiga())+"g\""
    if (pathBed == "")
      """
      ${params.path_gatk} --java-options ${java_options} HaplotypeCaller -I ${sample}.bam -O ${sample}_gatkHC_nofilter.vcf -R ${pathFasta} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --native-pair-hmm-threads ${task.cpus} --verbosity ERROR --max-reads-per-alignment-start 0
      """
    else
      """
      ${params.path_gatk} --java-options ${java_options} HaplotypeCaller -I ${sample}.bam -O ${sample}_gatkHC.vcf -R ${pathFasta} -L ${pathBed} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --native-pair-hmm-threads ${task.cpus} --verbosity ERROR --max-reads-per-alignment-start 0
      """
}


/* Mutect2 (only for mitochondrial mode)*/
process Mutect2Caller {
    label 'cpus_max' // allow group process with same computing reuqirements

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam) from ch_SampleToMutect2
      file bamvalidated from ch_BamValidated

    output:
      file("${sample}_mutect2.vcf") into ch_vcfMutect2

    when: (test==false) //mitoMode==true

    script:
    java_options = "\"-Xms"+(task.memory.toGiga() / 2).trunc()+"g -Xmx"+(task.memory.toGiga())+"g\""
    if (pathBed == "")
      """
      ${params.path_gatk} --java-options ${java_options} Mutect2 -I ${sample}.bam -O ${sample}_mutect2.vcf -R ${pathFasta} -min-AF ${params.min_af} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --initial-tumor-lod 0 --tumor-lod-to-emit 0 --af-of-alleles-not-in-resource 4e-3 --pruning-lod-threshold -4 --independent-mates
        sed -i s/"INFO=<ID=AS_SB_TABLE,Number=1"/"INFO=<ID=AS_SB_TABLE,Number=."/ ${sample}_mutect2.vcf
      """
    else
      """
      ${params.path_gatk} --java-options ${java_options} Mutect2 -I ${sample}.bam -O ${sample}_mutect2.vcf -R ${pathFasta} -L ${pathBed} -min-AF ${params.min_af} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --initial-tumor-lod 0 --tumor-lod-to-emit 0 --af-of-alleles-not-in-resource 4e-3 --pruning-lod-threshold -4 --independent-mates
        sed -i s/"INFO=<ID=AS_SB_TABLE,Number=1"/"INFO=<ID=AS_SB_TABLE,Number=."/ ${sample}_mutect2.vcf
      """
}



/* Deepvariant */
process Deepvariant {
    label 'cpus_max' // allow group process with same computing reuqirements

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam) from ch_SampleToDeepvariant
      file bamvalidated from ch_BamValidated

    output:
      file("${sample}_deepvariant.vcf") into ch_vcfDeepvariant

    when: (test==false) //platform=="IonTorrent"

    // DeepVariant require FASTA reference + BED target + BAM in the same folder => hard link files
    script:
      if (pathBed == "")
        """
        mkdir deepvariant_hl
        ln -f ${pathFasta} deepvariant_hl/${params.genome}.fasta
        ln -f ${pathFasta}.fai deepvariant_hl/${params.genome}.fasta.fai
        ln -f ${params.path_bam}/${sample}.bam deepvariant_hl/${sample}.bam
        ln -f ${params.path_bam}/${sample}.bam.bai deepvariant_hl/${sample}.bam.bai
        docker run -v \"\$PWD/deepvariant_hl\":\"/input\" -v \"\$PWD\":\"/output\" gcr.io/deepvariant-docker/deepvariant:\"${params.version_deepvariant}\" /opt/deepvariant/bin/run_deepvariant --model_type=WES --num_shards=${task.cpus} --ref=\"/input/${params.genome}.fasta\" --reads=\"/input/${sample}.bam\" --output_vcf=\"/output/${sample}_deepvariant.vcf.gz\" --novcf_stats_report
        rm -rf deepvariant_hl
        gzip -d ${sample}_deepvariant.vcf.gz ; rm -f /${sample}_deepvariant.vcf.gz.tbi
        """
      else
        """
        mkdir deepvariant_hl
        ln -f ${pathFasta} deepvariant_hl/${params.genome}.fasta
        ln -f ${pathFasta}.fai deepvariant_hl/${params.genome}.fasta.fai
        ln -f ${params.path_bam}/${sample}.bam deepvariant_hl/${sample}.bam
        ln -f ${params.path_bam}/${sample}.bam.bai deepvariant_hl/${sample}.bam.bai
        ln -f ${pathBed} deepvariant_hl/${targetName}
        docker run -v \"\$PWD/deepvariant_hl\":\"/input\" -v \"\$PWD\":\"/output\" gcr.io/deepvariant-docker/deepvariant:\"${params.version_deepvariant}\" /opt/deepvariant/bin/run_deepvariant --model_type=WES --num_shards=${task.cpus} --ref=\"/input/${params.genome}.fasta\" --reads=\"/input/${sample}.bam\" --regions=\"/input/${targetName}\" --output_vcf=\"/output/${sample}_deepvariant.vcf.gz\"
        rm -rf deepvariant_hl
        gzip -d ${sample}_deepvariant.vcf.gz ; rm -f ${sample}_deepvariant.vcf.gz.tbi
        """
}


/* Torrent Variant Caller (only for IonTorrent platform) */
process TorrentVariantCaller {
    label 'cpus_max' // allow group process with same computing reuqirements

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam) from ch_SampleToTVC
      file bamvalidated from ch_BamValidated

    output:
      file("${sample}_tvc.vcf") into ch_vcfTVC

    when: (test==true) //platform=="IonTorrent"

    script:
      if (pathBed == "")
        """
        ${params.path_tvc} -r ${pathFasta} -b ${sample}.bam -O \$PWD --output-vcf ${sample}_tvc.vcf --assembly-vcf ${sample}_tvc.vcf -M ${params.min_mapq} --num-threads ${task.cpus} --gen-min-alt-allele-freq ${params.min_af} --gen-min-indel-alt-allele-freq ${params.min_af} --snp-min-allele-freq ${params.min_af} --mnp-min-allele-freq ${params.min_af} --indel-min-allele-freq ${params.min_af} --gen-min-coverage ${params.min_varcov} --snp-min-var-coverage ${params.min_varcov} --mnp-min-var-coverage ${params.min_varcov} --indel-min-var-coverage ${params.min_varcov} --snp-min-variant-score ${params.min_varscore} --mnp-min-variant-score ${params.min_varscore} --indel-min-variant-score ${params.min_varscore} --snp-min-coverage ${params.min_cov} --mnp-min-coverage ${params.min_cov} --indel-min-coverage ${params.min_cov} --snp-strand-bias ${params.max_sb} --mnp-strand-bias ${params.max_sb} --indel-strand-bias ${params.max_sb} --allow-complex on --error-motifs-dir ${workflow.projectDir}/conf --error-motifs TVC_motifset.txt
        perl $workflow.projectDir/scripts/vcf-concat.pl -p ${sample}_tvc_filtered.vcf ${sample}_tvc.vcf > ${sample}_tvc_cat.vcf ; mv ${sample}_tvc_cat.vcf ${sample}_tvc.vcf
        """
      else
        """
        ${params.path_tvc} -r ${pathFasta} -b ${sample}.bam -t ${pathBed} -O \$PWD --output-vcf ${sample}_tvc.vcf --assembly-vcf ${sample}_tvc.vcf -M ${params.min_mapq} --num-threads ${task.cpus} --gen-min-alt-allele-freq ${params.min_af} --gen-min-indel-alt-allele-freq ${params.min_af} --snp-min-allele-freq ${params.min_af} --mnp-min-allele-freq ${params.min_af} --indel-min-allele-freq ${params.min_af} --gen-min-coverage ${params.min_varcov} --snp-min-var-coverage ${params.min_varcov} --mnp-min-var-coverage ${params.min_varcov} --indel-min-var-coverage ${params.min_varcov} --snp-min-variant-score ${params.min_varscore} --mnp-min-variant-score ${params.min_varscore} --indel-min-variant-score ${params.min_varscore} --snp-min-coverage ${params.min_cov} --mnp-min-coverage ${params.min_cov} --indel-min-coverage ${params.min_cov} --snp-strand-bias ${params.max_sb} --mnp-strand-bias ${params.max_sb} --indel-strand-bias ${params.max_sb} --allow-complex on --error-motifs-dir ${workflow.projectDir}/conf --error-motifs TVC_motifset.txt
        perl $workflow.projectDir/scripts/vcf-concat.pl -p ${sample}_tvc_filtered.vcf ${sample}_tvc.vcf > ${sample}_tvc_cat.vcf ; mv ${sample}_tvc_cat.vcf ${sample}_tvc.vcf
        """
}


/* Strelka2 (only for Illumina platform) */
process Strelka2 {
    label 'cpus_max' // allow group process with same computing reuqirements

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam) from ch_SampleToStrelka2
      file bamvalidated from ch_BamValidated

    output:
      file("${sample}_strelka.vcf") into ch_vcfStrelka2

    when: (test==false) //platform=="Illumina"

    script:
      // Strelka require a sort, tabix and gzip BED file => pathBedgz
      if (pathBed == "")
        """
        python2.7 ${params.path_strelka} --bam ${sample}.bam --referenceFasta ${pathFasta} -targeted --runDir ${sample}_strelka
        python2.7 -c \"import pickle ; data=pickle.load(open(\\\""+dicoNiourk["tmp_dir"]+"/"+bc+"_strelka/runWorkflow.py.config.pickle\\\",\\\"rb\\\")) data['StrelkaGermline']['minMapq']=${params.min_mapq} ; pickle.dump(data,open(\\\""+dicoNiourk["tmp_dir"]+"/"+bc+"_strelka/runWorkflow.py.config.pickle\\\", \\\"w\\\"))\"
        sed -i '1i # -*- coding: utf-8 -*-' ${sample}_strelka/runWorkflow.py
        python2.7 ${sample}_strelka/runWorkflow.py -m local -j ${task.cpus} -g ${task.memory.toGiga()}
        mv ${sample}_strelka/results/variants/variants.vcf.gz ${sample}_strelka.vcf.gz
        gzip -d ${sample}_strelka.vcf.gz
        """
      else
        """
        python2.7 ${params.path_strelka} --bam ${sample}.bam --referenceFasta ${pathFasta} --targeted --runDir ${sample}_strelka --callRegions ${workflow.projectDir}/${pathBedgz}
        python2.7 -c \"import pickle ; data=pickle.load(open(\\\"${sample}_strelka/runWorkflow.py.config.pickle\\\",\\\"rb\\\")) ; data['StrelkaGermline']['minMapq']=${params.min_mapq} ; pickle.dump(data,open(\\\"${sample}_strelka/runWorkflow.py.config.pickle\\\", \\\"w\\\"))\"
        sed -i '1i # -*- coding: utf-8 -*-' ${sample}_strelka/runWorkflow.py
        python2.7 ${sample}_strelka/runWorkflow.py -m local -j ${task.cpus} -g ${task.memory.toGiga()}
        mv ${sample}_strelka/results/variants/variants.vcf.gz ${sample}_strelka.vcf.gz
        gzip -d ${sample}_strelka.vcf.gz
        """
}

/*      
**********************************************************************
                      VARIANTS PROCESSING                    
**********************************************************************
*/
/* Merge VCFs files (format, normalize, decompose, merge) */
process MergeSampleCallsVCFs {
    label 'cpus_8' // allow group process with same computing reuqirements

    tag "${sample}" // tag in report tasks

    publishDir pathVCF, mode: 'copy', overwrite: true
    
    input:
      set sample, file(bam) from ch_SampleVCF
      file vcf from ch_vcfGATKhc.mix(ch_vcfMutect2,ch_vcfDeepvariant,ch_vcfTVC,ch_vcfStrelka2).collect()

    output:
      set sample, file("${sample}_NkMerged.vcf") into ch_MergedVCF

    script:
    """
    python3 $workflow.projectDir/scripts/Nk_mergeVCFs.py ${params.path_gatk} ${params.path_vt} ${params.path_vcfvalidator} ${pathVCFraw} ${pathFasta} ${params.min_cov} ${sample} ${vcf}
    """
}


/* VEP annotation */
process VEPannotateVCFs {
    label 'cpus_8' // allow group process with same computing reuqirements

    tag "${sample}" // tag in report tasks

    publishDir pathVCF, mode: 'copy', overwrite: true
    
    input:
      set sample, file(mergedVcf) from ch_MergedVCF

    output:
      file("${sample}_NkMergedVEP.vcf") into ch_VEPVCF

    script:
    """
    ${params.path_vep} --offline --cache --dir ${params.path_vep_cache} --refseq --species homo_sapiens --assembly ${params.genome} --input_file ${mergedVcf} --output_file ${sample}_NkMerged_vep.vcf --force_overwrite --stats_file ${sample}_NkMerged_vep_summary.html --fork 8 --everything --hgvs --check_ref --fasta ${pathFasta} --vcf
    """
}

