#!/usr/bin/env nextflow


/*
**********************************************************************
                             HELP & USAGE                             
**********************************************************************
*/
// Colors
c_reset  = params.monochrome ? '' : "\033[0m";
c_dim    = params.monochrome ? '' : "\033[2m";
c_header = params.monochrome ? '' : "\033[38;2;93;191;230m";
c_subtitle = params.monochrome ? '' : "\033[38;2;230;132;93m";
c_light1 = params.monochrome ? '' : "\033[38;2;252;146;86m";
c_light2 = params.monochrome ? '' : "\033[38;2;200;200;200m";
c_grey = params.monochrome ? '' : "\033[38;2;150;150;150m";
c_red = params.monochrome ? '' : "\033[38;2;212;64;89m";
c_green = params.monochrome ? '' : "\033[38;2;55;200;113m";

def helpMessage(strError) {
    log.info niourkHeader()
    println c_reset
    log.error"""
    $strError
    """
    println c_reset
    log.info"""
    🆄 🆂 🅰 🅶 🅴
    The typical command for running the pipeline is as follows:

      nextflow run niourk.nf --pathBam bamfolder --pathOut results --genome GRCh37 [arguments]
      
      Mandatory arguments:
        --path_bam                  Path to input BAMs folder
        --path_out                  Path to results folder
        --genome                    Reference Genome (GRCh37, GRCh38, rCRS, ...)
        --mito                      Mitochondrial mode (true|false)

      Sequencing arguments:
        --wgs                       Whole Genome Mode
        --path_bed                  Path to target BED file
        --path_genes                Path to genes list file (one gene per line)
        --platform                  Sequencing platform (iontorrent|illumina)

      Tools arguments:
        --path_elprep               Path to elPrep executable
        --path_elprep_files         Path to elPrep required files (elsites)
        --path_gatk                 Path to Genome Analysis Toolkit (GATK) executable
        --path_tvc                  Path to Torrent Variant Caller (TVC) executable
        --path_strelka              Path to Strelka executable
        --path_samtools             Path to samtools executable
        --path_mosdepth             Path to mosdepth executable
        --path_vt                   Path to vt executable
        --path_vep                  Path to VEP source directory
        --path_vep_cache            Path to VEP cache directory
        --path_vep_plugin_files     Path to VEP plugins files
        --path_vcfanno              Path to VcfAnno executable
        --path_vcfvalidator         Path to EBIvariation vcf-validator
        --version_deepvariant       Version tag for deepvariant docker container

      Depth arguments:
        --padding                   Padding size to consider

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
        --monochrome                Disable ansi colors
    """.stripIndent()
    println c_reset
}
// Show help message
if (params.help) exit 0, helpMessage("")

// Check input BAMs
if (!params.path_bam) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nParameter `--path_bam` is required")
if (!file(params.path_bam).exists()) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid directory for `--path_bam` (`"+params.path_bam+"`)")
bamFiles = file(params.path_bam).list()
bamCount = 0
for( def bam : bamFiles ) { if (file(bam).getExtension()=="bam") bamCount+=1 }
if (bamCount==0) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nAny BAM file found in `--path_bam` directory (`"+params.path_bam+"`)")

// Check output BAMs
if (!params.path_out) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nParameter `--path_out` is required")
if (!file(params.path_out).isDirectory()) file(params.path_out).mkdir()
pathDepth= file(params.path_out+"/depth")
if (!pathDepth.isDirectory()) pathDepth.mkdir()
pathVCF = file(params.path_out+"/vcf")
if (!pathVCF.isDirectory()) pathVCF.mkdir()
pathVCFraw = file(params.path_out+"/vcf/raw")
if (!pathVCFraw.isDirectory()) pathVCFraw.mkdir()
pathFile = file(params.path_out+"/files")
if (!pathFile.isDirectory()) pathFile.mkdir()

// Check reference
if (params.genome=="") exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nParameter `--genome` is required")
pathGenomeDir = "${workflow.projectDir}/reference/${params.genome}"
if (!file(pathGenomeDir).isDirectory()) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid genome assembly:`GRCh37` or `GRCh38`")
pathFasta = file(pathGenomeDir+"/"+params.genome+".fasta")
pathFastaFai = file(pathGenomeDir+"/"+params.genome+".fasta.fai")
pathFastaElprep = file(pathGenomeDir+"/"+params.genome+".elfasta")
if (!pathFasta.exists()) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nAny FASTA file found for specified `--genome` (`"+params.genome+"`)")
if (!pathFastaFai.exists()) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nAny FASTA index file found for specified `--genome` (`"+params.genome+"`)")
if (params.genome=="GRCh37") ponp2_ref = "hg37"
if (params.genome=="GRCh38") ponp2_ref = "hg38"

// GFF & gffutils DB
pathGff = file(pathGenomeDir+"/"+params.genome+".gff")
pathGffDB = file(pathGenomeDir+"/"+params.genome+".gff.db")
if (!pathGff.exists()) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nAny GFF file found for specified `--genome` (`"+params.genome+"`)")
if (!pathGffDB.exists()) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nAny GFF gffutils database file found for specified `--genome` (`"+params.genome+"`)")

// VEP Cache directory
if (!file(params.path_vep_cache).isDirectory())  exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid VEP cache directory for `--path_vep_cache` (`"+params.path_vep_cache+"`)")
if (!file(params.path_vep_cache+"/homo_sapiens_refseq").isDirectory())  exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nAny `homo_sapiens_refseq` directory found in VEP cache directory for `--path_vep_cache` (`"+params.path_vep_cache+"`)")
if (!file(params.path_vep_cache+"/homo_sapiens_refseq/100_"+params.genome).isDirectory())  exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nAny `100_"+params.genome+"` directory found in VEP cache subdirectory for `--path_vep_cache` (`"+params.path_vep_cache+"/homo_sapiens_refseq`)")

// Mitochondrial mode
mitoMode = params.mito.toBoolean()
if ((params.mito!=false) && (params.mito!=true)) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid mitochondrial mode `"+params.mito+"`: `true` or `false`")

// WGS mode
wgsMode = params.wgs.toBoolean()
if ((params.wgs!=false) && (params.wgs!=true)) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid wgs mode `"+params.mito+"`: `true` or `false`")

// Sequencing platform
if (params.platform=="") exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nParameter `--platform` is required")
if ((params.platform!="iontorrent") && (params.platform!="illumina")) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid sequencing platform:`iontorrent` or `illumina`")

// Target BED file
pathBed = params.path_bed
if (pathBed!="")
  {
  if ((pathBed!=null) && (!file(pathBed).exists())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid file for `--path_bed` (`"+params.path_bed+"`)")
  }
else { pathBed = "None" }

// Genes list file
pathGenes = params.path_genes
if (pathGenes!="")
  {
  if ((pathGenes!=null) && (!file(pathGenes).exists())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid file for `--path_genes` (`"+params.path_genes+"`)")
  }
else { pathGenes = "None" }

// Both Target BED and genes list file error
if ((pathBed!="None") && (pathGenes!="None")) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nBoth `--path_bed` and `--path_genes` is not allowed, please choose one")

// NiourK BED file
if (wgsMode==false)
  {
   pathNiourkBed = pathFile+"/NiourK.bed"
   pathNiourkBedgz = pathFile+"/NiourK.bed.gz"
  } 
else{
   pathNiourkBed = ""
   pathNiourkBedgz = ""
  } 
// Check tools path
missing_tools = ""
if ((params.path_elprep=="") || (!file(params.path_elprep).exists())) missing_tools+="\n    `--path_elprep`"
if ((params.path_elprep_files=="") || (!file(params.path_elprep_files).isDirectory())) missing_tools+="\n   `--path_elprep_files`"
if ((params.path_gatk=="") || (!file(params.path_gatk).exists())) missing_tools+="\n    `--path_gatk`"
if ((params.path_tvc=="") || (!file(params.path_tvc).exists())) missing_tools+="\n   `--path_tvc`"
if ((params.path_strelka=="") || (!file(params.path_strelka).exists())) missing_tools+="\n   `--path_strelka`"
if ((params.path_samtools=="") || (!file(params.path_samtools).exists())) missing_tools+="\n   `--path_samtools`"
if ((params.path_mosdepth=="") || (!file(params.path_mosdepth).exists())) missing_tools+="\n   `--path_mosdepth`"
if ((params.path_vt=="") || (!file(params.path_vt).exists())) missing_tools+="\n   `--path_vt`"
if ((params.path_vep=="") || (!file(params.path_vep).exists())) missing_tools+="\n   `--path_vep`"
if ((params.path_vep_cache=="") || (!file(params.path_vep_cache).isDirectory())) missing_tools+="\n   `--path_vep_cache`"
if ((params.path_vep_plugin_files=="") || (!file(params.path_vep_plugin_files).isDirectory())) missing_tools+="\n   `--path_vep_plugin_files`"
if ((params.path_vcfanno=="") || (!file(params.path_vcfanno).exists())) missing_tools+="\n   `--path_vcfanno`"
if ((params.path_vcfvalidator=="") || (!file(params.path_vcfvalidator).exists())) missing_tools+="\n   `--path_vcfvalidator`"
if (missing_tools!="") exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid path for params: "+missing_tools)
if ((params.version_deepvariant=="") || (!params.version_deepvariant)) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid version for `--version_deepvariant` (`"+params.version_deepvariant+"`)")

// Depth
if ((params.padding=="") || (!params.padding.toString().isNumber())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid integer for `--padding` (`"+params.padding.toString()+"`)")

// Calling
if ((params.min_baseq=="") || (!params.min_baseq.toString().isNumber())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid integer for `--min_baseq` (`"+params.min_baseq.toString()+"`)")
if ((params.min_mapq=="") || (!params.min_mapq.toString().isNumber())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid integer for `--min_mapq` (`"+params.min_mapq.toString()+"`)")
if ((params.min_af=="") || (!params.min_af.toString().isFloat())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid float for `--min_af` (`"+params.min_af.toString()+"`)")
if ((params.min_cov=="") || (!params.min_cov.toString().isNumber())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid integer for `--min_cov` (`"+params.min_cov.toString()+"`)")
if ((params.min_varcov=="") || (!params.min_varcov.toString().isNumber())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid integer for `--min_varcov` (`"+params.min_varcov.toString()+"`)")
if ((params.min_varscore=="") || (!params.min_varscore.toString().isNumber())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid integer for `--min_varscore` (`"+params.min_varscore.toString()+"`)")
if ((params.max_sb=="") || (!params.max_sb.toString().isFloat())) exit 1, helpMessage("\n🅴 🆁 🆁 🅾 🆁\nPlease specify a valid float for `--max_sb` (`"+params.max_sb.toString()+"`)")

// Channels
Channel
    .fromFilePairs(params.path_bam+'/*.{bam,bai}') { file -> file.name.replaceAll(/.bam|.bai$/,'') }
    .set { ch_Sample }

// split bam channel into multiples
ch_Sample.into { ch_SampleToValidate ; ch_SampleToPreprocess }



/*
**********************************************************************
                           DISPLAY & SUMMARY                          
**********************************************************************
*/
/* Summary */
log.info niourkHeader() ; println c_reset
log.info digestHeader() ; println c_reset
println "🅴 🆇 🅴 🅲 🆄 🆃 🅸 🅾 🅽"

/* ABOUT header */
def niourkHeader() {""
    if (!params.monochrome) println c_header
    Date start = new Date();
    startFormat = start.toString(); 
    return """
🇳​​​​​   🇮​​​​​   🇴​​​​​   🇺​​​​​   🇷​​​​​   🇰​​

🅰 🅱 🅾 🆄 🆃
Workflow   : NGS for clinical diagnostics
Version    : v${workflow.manifest.version}
Copyright  : KAVE - CHU Angers
Support    : https://github.com/dooguypapua/NiourK_nf
    """.stripIndent()
}

/*  DIGEST header */
def digestHeader() {""
    if (!params.monochrome) print c_subtitle
    if (params.platform=="IonTorrent")
    {
      if (mitoMode==true) callers = "TVC, Mutect2, Deepvariant"
      else callers = "TVC, HaplotypeCaller, Deepvariant"
    }
    else callers = "Strelka2, HaplotypeCaller, Deepvariant"
    return """
🅳 🅸 🅶 🅴 🆂 🆃
Date       : ${startFormat}
Hostname   : ${InetAddress.localHost.canonicalHostName}
User       : ${workflow.userName}
Cpus       : ${params.max_cpus}
Memory     : ${params.max_memory}
NfDir      : ${workflow.launchDir}
WorkDir    : ${workflow.workDir}
Report     : ${workflow.sessionId}
Input      : ${params.path_bam}
Output     : ${params.path_out}
Genome     : ${params.genome}
BAMs       : ${bamCount} file(s) found
Target     : ${pathBed}
Platform   : ${params.platform}
MitoMode   : ${mitoMode}
Callers    : ${callers}
VEPplugins : Blosum62, Conservation, dbscSNV, CADD, REVEL, Mastermind, SpliceAI
    """.stripIndent()
}

/*
**********************************************************************
                       CHECKING REQUIRED FILES                        
**********************************************************************
*/
/* Create NiourK bed file */
process MakeNiourkBed {
    cache 'deep'
    script:
    if (pathNiourkBedgz != "")
      """
      python3 $workflow.projectDir/scripts/Nk_niourkBed.py $pathBed $pathGenes $pathGff $pathFile ${params.padding}
      """
}

/* Parse software version numbers */
process GetSoftwareVersions {
    cache 'deep'
    script:
      """
      python3 $workflow.projectDir/scripts/Nk_checkTools.py gatk=${params.path_gatk} tvc=${params.path_tvc} \
strelka=${params.path_strelka} samtools=${params.path_samtools} vt=${params.path_vt} \
mosdepth=${params.path_mosdepth} vcfanno=${params.path_vcfanno} vcf-validator=${params.path_vcfvalidator} \
deepvariant=${params.version_deepvariant} vep=${params.path_vep}/vep vepplugin=${params.path_vep_plugin_files} > ${pathFile}/tools_version.json
      """
}


/* Validate & index input BAM files */
process ValidateBamFile {
    cache 'deep'
    label 'cpus_8' // allow group process with same computing requirements

    tag "${sample}" // tag in report tasks

    errorStrategy "terminate"

    input:
      set sample, file(bam) from ch_SampleToValidate

    output:
      file ("BamValidated.txt") into ch_BamValidated

    script:
      java_options = "\"-Xms"+(task.memory.toGiga() / 2).trunc()+"g -Xmx"+(task.memory.toGiga())+"g\""
      """
      ${params.path_gatk} --java-options ${java_options} ValidateSamFile -I ${sample}.bam >/dev/null 2>&1 || (>&2 echo "🅴 🆁 🆁 🅾 🆁 ValidateSamFile" && exit 1)

      python3 $workflow.projectDir/scripts/Nk_checkBam.py ${params.path_samtools} ${task.cpus} ${sample}.bam $workflow.projectDir/reference/${params.genome}/${params.genome}.dict

      echo "SUCCESS" > BamValidated.txt
      """
}



/*      
**********************************************************************
                           BAM MANIPULATION                           
**********************************************************************
*/
/* BAM preprocessing with elprep (BaseRecalibrator+ApplyBQSR) 
   For iontorrent = simplify BAM header before elprep and restore header after
*/
process BamPreprocessing {
    cache 'deep'
    label 'cpus_max' // allow group process with same computing requirements
    label 'memory_max'
    publishDir pathFile, mode: 'copy', overwrite: true

    input: 
      set sample, file(bam) from ch_SampleToPreprocess
      file bamvalidated from ch_BamValidated

    output:
      set sample, file("${sample}_elprep.bam"), file("${sample}_elprep.bam.bai") into ch_BamPreprocessed
    
    script:    
      if (params.platform=="iontorrent")
      """
      samtools view -H ${sample}.bam > ${sample}_header.sam
      samtools view -@ 8 -h ${sample}.bam | grep -Pv "@[CO|PV]" | ${params.path_elprep} sfm /dev/stdin /dev/stdout --filter-unmapped-reads --filter-mapping-quality ${params.min_mapq} --mark-duplicates --mark-optical-duplicates ${sample}_mark-optical-duplicates.txt --optical-duplicates-pixel-distance 2500 --remove-duplicates --sorting-order coordinate --clean-sam --bqsr ${sample}_bqsr_recall.txt --bqsr-reference ${pathFastaElprep} --known-sites ${params.path_elprep_files}/1000G-known-indels_${params.genome}.elsites,${params.path_elprep_files}/1000G-omni2.5_${params.genome}.elsites,${params.path_elprep_files}/1000G-phase1-snps-highconfidence_${params.genome}.elsites,${params.path_elprep_files}/dbsnp138_${params.genome}.elsites,${params.path_elprep_files}/Millsand1000G-goldstandard-indels_${params.genome}.elsites --nr-of-threads 8 --log-path ${sample}_elprepLog --intermediate-files-output-prefix ${sample} --intermediate-files-output-type sam --tmp-path /tmp | samtools view -@ 8 -b | samtools reheader ${sample}_header.sam - > ${sample}_elprep.bam
      samtools index ${sample}_elprep.bam
      cp ${sample}_mark-optical-duplicates.txt ${pathFile}
      cp ${sample}_bqsr_recall.txt ${pathFile}
      """
      else
      """
      ${params.path_elprep} sfm ${sample}.bam ${sample}_elprep.bam --filter-unmapped-reads --filter-mapping-quality ${params.min_mapq} --mark-duplicates --mark-optical-duplicates ${sample}_mark-optical-duplicates.txt --optical-duplicates-pixel-distance 2500 --remove-duplicates --sorting-order coordinate --clean-sam --bqsr ${sample}_bqsr_recall.txt --bqsr-reference ${pathFastaElprep} --known-sites ${params.path_elprep_files}/1000G-known-indels_${params.genome}.elsites,${params.path_elprep_files}/1000G-omni2.5_${params.genome}.elsites,${params.path_elprep_files}/1000G-phase1-snps-highconfidence_${params.genome}.elsites,${params.path_elprep_files}/dbsnp138_${params.genome}.elsites,${params.path_elprep_files}/Millsand1000G-goldstandard-indels_${params.genome}.elsites --nr-of-threads 8 --log-path ${sample}_elprepLog --intermediate-files-output-prefix ${sample} --intermediate-files-output-type sam --tmp-path /tmp
      samtools index ${sample}_elprep.bam
      cp ${sample}_mark-optical-duplicates.txt ${pathFile}
      cp ${sample}_bqsr_recall.txt ${pathFile}
      """
}

// split bam channel into multiples
ch_BamPreprocessed.into { ch_SampleToDepth ; ch_SampleToGATKhc ; ch_SampleToMutect2 ; ch_SampleToTVC ; ch_SampleToDeepvariant ; ch_SampleToStrelka2 ; ch_SampleVCF }


/* Compute BAM depth with mosdepth */
process BamDepth {
    cache 'deep'
    label 'cpus_8' // allow group process with same computing requirements

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam), file(bai) from ch_SampleToDepth

    output:
      file("${sample}.mosdepth*") into ch_bamDepth

    script:
      /* --thresholds 10,20,50 */
      """
      ${params.path_mosdepth} -t ${task.cpus} -b ${pathNiourkBedgz} -Q ${params.min_mapq} --fast-mode --use-median ${sample} ${bam}
      cp ${sample}.mosdepth* ${pathDepth}
      cp ${sample}*.bed.gz* ${pathDepth}
      """
}



/*      
**********************************************************************
                               CALLING                                
**********************************************************************
// All VCFs are created in their work dir before validation and copy to pathVCFraw
*/
/* Deepvariant */
process Deepvariant {
    cache 'deep'
    label 'cpus_max' // allow group process with same computing requirements
    label 'memory_max'

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam), file(bai) from ch_SampleToDeepvariant

    output:
      file("${sample}_deepvariant.vcf") into ch_vcfDeepvariant

    // DeepVariant require FASTA reference + BED target + BAM in the same folder => hard link files
    script:
      if (pathNiourkBedgz == "")
        """
        mkdir deepvariant_hl
        ln -f ${pathFasta} deepvariant_hl/${params.genome}.fasta
        ln -f ${pathFasta}.fai deepvariant_hl/${params.genome}.fasta.fai
        ln -f ${pathFile}/${sample}_elprep.bam deepvariant_hl/${sample}_elprep.bam
        ln -f ${pathFile}/${sample}_elprep.bam.bai  deepvariant_hl/${sample}_elprep.bam.bai
        docker run -v \"\$PWD/deepvariant_hl\":\"/input\" -v \"\$PWD\":\"/output\" gcr.io/deepvariant-docker/deepvariant:\"${params.version_deepvariant}\" /opt/deepvariant/bin/run_deepvariant --model_type=WES --num_shards=${task.cpus} --ref=\"/input/${params.genome}.fasta\" --reads=\"/input/${sample}_elprep.bam\" --output_vcf=\"/output/${sample}_deepvariant.vcf.gz\" --novcf_stats_report
        rm -rf deepvariant_hl
        gzip -d ${sample}_deepvariant.vcf.gz ; rm -f /${sample}_deepvariant.vcf.gz.tbi
        """
      else
        """
        mkdir deepvariant_hl
        ln -f ${pathFasta} deepvariant_hl/${params.genome}.fasta
        ln -f ${pathFasta}.fai deepvariant_hl/${params.genome}.fasta.fai
        ln -f ${pathFile}/${sample}_elprep.bam deepvariant_hl/${sample}_elprep.bam
        ln -f ${pathFile}/${sample}_elprep.bam.bai  deepvariant_hl/${sample}_elprep.bam.bai
        gzip -d -k -c ${pathNiourkBedgz} | cut -f 1,2,3 > deepvariant_hl/NiourK.bed
        docker run -v \"\$PWD/deepvariant_hl\":\"/input\" -v \"\$PWD\":\"/output\" gcr.io/deepvariant-docker/deepvariant:\"${params.version_deepvariant}\" /opt/deepvariant/bin/run_deepvariant --model_type=WES --num_shards=${task.cpus} --ref=\"/input/${params.genome}.fasta\" --reads=\"/input/${sample}_elprep.bam\" --regions=\"/input/NiourK.bed\" --output_vcf=\"/output/${sample}_deepvariant.vcf.gz\"
        rm -rf deepvariant_hl
        gzip -d ${sample}_deepvariant.vcf.gz ; rm -f ${sample}_deepvariant.vcf.gz.tbi
        """
}


/* GATK HaplotypeCaller */
process HaplotypeCaller {
    cache 'deep'
    label 'cpus_max' // allow group process with same computing requirements
    label 'memory_max'

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam), file(bai) from ch_SampleToGATKhc

    output:
      file("${sample}_gatkHC.vcf") into ch_vcfGATKhc

    when: (mitoMode==false)

    script:
    java_options = "\"-Xms"+(task.memory.toGiga() / 2).trunc()+"g -Xmx"+(task.memory.toGiga())+"g\""
    if (pathNiourkBedgz == "")
      """
      ${params.path_gatk} --java-options ${java_options} HaplotypeCaller -I ${bam} -O ${sample}_gatkHC.vcf -R ${pathFasta} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --native-pair-hmm-threads ${task.cpus} --verbosity ERROR --max-reads-per-alignment-start 0
      """
    else
      """
      ${params.path_gatk} --java-options ${java_options} HaplotypeCaller -I ${bam} -O ${sample}_gatkHC.vcf -R ${pathFasta} -L ${pathNiourkBedgz} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --native-pair-hmm-threads ${task.cpus} --verbosity ERROR --max-reads-per-alignment-start 0
      """
}


/* Mutect2 (only for mitochondrial mode)*/
process Mutect2Caller {
    cache 'deep'
    label 'cpus_max' // allow group process with same computing requirements
    label 'memory_max'

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam), file(bai) from ch_SampleToMutect2

    output:
      file("${sample}_mutect2.vcf") into ch_vcfMutect2

    when: (mitoMode==true)

    script:
    java_options = "\"-Xms"+(task.memory.toGiga() / 2).trunc()+"g -Xmx"+(task.memory.toGiga())+"g\""
    if (pathNiourkBedgz == "")
      """
      ${params.path_gatk} --java-options ${java_options} Mutect2 -I ${bam} -O ${sample}_mutect2.vcf -R ${pathFasta} -min-AF ${params.min_af} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --initial-tumor-lod 0 --tumor-lod-to-emit 0 --af-of-alleles-not-in-resource 4e-3 --pruning-lod-threshold -4 --independent-mates
        sed -i s/"INFO=<ID=AS_SB_TABLE,Number=1"/"INFO=<ID=AS_SB_TABLE,Number=."/ ${sample}_mutect2.vcf
      """
    else
      """
      ${params.path_gatk} --java-options ${java_options} Mutect2 -I ${bam} -O ${sample}_mutect2.vcf -R ${pathFasta} -L ${pathNiourkBedgz} -min-AF ${params.min_af} -mbq ${params.min_baseq} --minimum-mapping-quality ${params.min_mapq} --initial-tumor-lod 0 --tumor-lod-to-emit 0 --af-of-alleles-not-in-resource 4e-3 --pruning-lod-threshold -4 --independent-mates
        sed -i s/"INFO=<ID=AS_SB_TABLE,Number=1"/"INFO=<ID=AS_SB_TABLE,Number=."/ ${sample}_mutect2.vcf
      """
}


/* Torrent Variant Caller (only for IonTorrent platform) */
process TorrentVariantCaller {
    cache 'deep'
    label 'cpus_max' // allow group process with same computing requirements
    label 'memory_max'

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam), file(bai) from ch_SampleToTVC

    output:
      file("${sample}_tvc.vcf") into ch_vcfTVC

    when: (params.platform=="iontorrent")

    script:
      if (pathNiourkBedgz == "")
        """
        ${params.path_tvc} -r ${pathFasta} -b ${bam} -O \$PWD --output-vcf ${sample}_tvc.vcf --assembly-vcf ${sample}_tvc.vcf -M ${params.min_mapq} --num-threads ${task.cpus} --gen-min-alt-allele-freq ${params.min_af} --gen-min-indel-alt-allele-freq ${params.min_af} --snp-min-allele-freq ${params.min_af} --mnp-min-allele-freq ${params.min_af} --indel-min-allele-freq ${params.min_af} --gen-min-coverage ${params.min_varcov} --snp-min-var-coverage ${params.min_varcov} --mnp-min-var-coverage ${params.min_varcov} --indel-min-var-coverage ${params.min_varcov} --snp-min-variant-score ${params.min_varscore} --mnp-min-variant-score ${params.min_varscore} --indel-min-variant-score ${params.min_varscore} --snp-min-coverage ${params.min_cov} --mnp-min-coverage ${params.min_cov} --indel-min-coverage ${params.min_cov} --snp-strand-bias ${params.max_sb} --mnp-strand-bias ${params.max_sb} --indel-strand-bias ${params.max_sb} --allow-complex on --error-motifs-dir ${workflow.projectDir}/conf --error-motifs TVC_motifset.txt
        perl $workflow.projectDir/scripts/vcf-concat.pl -p ${sample}_tvc_filtered.vcf ${sample}_tvc.vcf > ${sample}_tvc_cat.vcf ; mv ${sample}_tvc_cat.vcf ${sample}_tvc.vcf
        """
      else
        """
        ${params.path_tvc} -r ${pathFasta} -b ${bam} -t ${pathNiourkBed} -O \$PWD --output-vcf ${sample}_tvc.vcf --assembly-vcf ${sample}_tvc.vcf -M ${params.min_mapq} --num-threads ${task.cpus} --gen-min-alt-allele-freq ${params.min_af} --gen-min-indel-alt-allele-freq ${params.min_af} --snp-min-allele-freq ${params.min_af} --mnp-min-allele-freq ${params.min_af} --indel-min-allele-freq ${params.min_af} --gen-min-coverage ${params.min_varcov} --snp-min-var-coverage ${params.min_varcov} --mnp-min-var-coverage ${params.min_varcov} --indel-min-var-coverage ${params.min_varcov} --snp-min-variant-score ${params.min_varscore} --mnp-min-variant-score ${params.min_varscore} --indel-min-variant-score ${params.min_varscore} --snp-min-coverage ${params.min_cov} --mnp-min-coverage ${params.min_cov} --indel-min-coverage ${params.min_cov} --snp-strand-bias ${params.max_sb} --mnp-strand-bias ${params.max_sb} --indel-strand-bias ${params.max_sb} --allow-complex on --error-motifs-dir ${workflow.projectDir}/conf --error-motifs TVC_motifset.txt
        perl $workflow.projectDir/scripts/vcf-concat.pl -p ${sample}_tvc_filtered.vcf ${sample}_tvc.vcf > ${sample}_tvc_cat.vcf ; mv ${sample}_tvc_cat.vcf ${sample}_tvc.vcf
        """
}


/* Strelka2 (only for Illumina platform) */
process Strelka2 {
    cache 'deep'
    label 'cpus_max' // allow group process with same computing requirements
    label 'memory_max'

    tag "${sample}" // tag in report tasks

    input:
      set sample, file(bam), file(bai) from ch_SampleToStrelka2

    output:
      file("${sample}_strelka.vcf") into ch_vcfStrelka2

    when: (params.platform=="Illumina")

    script:
      // Strelka require a sort, tabix and gzip BED file => pathNiourkBedgzgz
      if (pathNiourkBedgz == "")
        """
        python2.7 ${params.path_strelka} --bam ${bam} --referenceFasta ${pathFasta} -targeted --runDir ${sample}_strelka
        python2.7 -c \"import pickle ; data=pickle.load(open(\\\""+dicoNiourk["tmp_dir"]+"/"+bc+"_strelka/runWorkflow.py.config.pickle\\\",\\\"rb\\\")) data['StrelkaGermline']['minMapq']=${params.min_mapq} ; pickle.dump(data,open(\\\""+dicoNiourk["tmp_dir"]+"/"+bc+"_strelka/runWorkflow.py.config.pickle\\\", \\\"w\\\"))\"
        sed -i '1i # -*- coding: utf-8 -*-' ${sample}_strelka/runWorkflow.py
        python2.7 ${sample}_strelka/runWorkflow.py -m local -j ${task.cpus} -g ${task.memory.toGiga()}
        mv ${sample}_strelka/results/variants/variants.vcf.gz ${sample}_strelka.vcf.gz
        gzip -d ${sample}_strelka.vcf.gz
        """
      else
        """
        python2.7 ${params.path_strelka} --bam ${bam} --referenceFasta ${pathFasta} --targeted --runDir ${sample}_strelka --callRegions ${workflow.projectDir}/${pathNiourkBedgzgz}
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
    cache 'deep'
    label 'cpus_8' // allow group process with same computing requirements
    label 'memory_max'

    tag "${sample}" // tag in report tasks

    publishDir pathVCF, mode: 'copy', overwrite: true
    
    input:
      set sample, file(bam) from ch_SampleVCF
      file vcf from ch_vcfGATKhc.mix(ch_vcfMutect2,ch_vcfDeepvariant,ch_vcfTVC,ch_vcfStrelka2).collect()

    output:
      set sample, file("${sample}_Nk.vcf.gz"), file("${sample}_Nk.vcf.gz.tbi") into ch_MergedVCF

    script:
    """
    python3 $workflow.projectDir/scripts/Nk_mergeVCFs.py ${params.path_gatk} ${params.path_vt} ${params.path_vcfvalidator} ${pathVCFraw} ${pathFasta} ${params.min_cov} ${sample} ${vcf}
    """
}


/* VEP annotation */
process VEPannotateVCFs {
    cache 'deep'
    label 'cpus_max' // allow group process with same computing requirements
    label 'memory_max'

    tag "${sample}" // tag in report tasks

    publishDir pathVCF, mode: 'copy', overwrite: true
    
    input:
      set sample, file(mergedVcf), file(mergedVcftbi) from ch_MergedVCF

    output:
      set sample, file("${sample}_Nk_vep.vcf.gz"),file("${sample}_Nk_vep.vcf.gz.tbi") into ch_VEPVCF

    script:
    """
    ${params.path_vep}/vep --offline --quiet --force_overwrite --cache --refseq --dir_cache ${params.path_vep_cache} --dir_plugins ${params.path_vep} \
--species homo_sapiens --assembly ${params.genome} --use_transcript_ref --check_ref --fasta ${pathFasta} \
--input_file ${mergedVcf} --output_file ${sample}_Nk_vep.vcf.gz  \
--vcf --compress_output bgzip --fork ${task.cpus} \
--sift p --polyphen p --af_gnomad --max_af --hgvs --variant_class --symbol --gene_phenotype --domains --numbers --biotype --pubmed \
--plugin Blosum62 \
--plugin Conservation,${params.path_vep_plugin_files}/GERP_${params.genome}.bw \
--plugin CADD,${params.path_vep_plugin_files}/CADD_${params.genome}.tsv.gz \
--plugin REVEL,${params.path_vep_plugin_files}/revel.tsv.gz \
--plugin Mastermind,${params.path_vep_plugin_files}/mastermind_${params.genome}.vcf.gz,1 \
--plugin dbscSNV,${params.path_vep_plugin_files}/dbscSNV_${params.genome}.txt.gz \
--plugin SpliceAI,snv=${params.path_vep_plugin_files}/spliceai-snv_${params.genome}.vcf.gz,indel=${params.path_vep_plugin_files}/spliceai-indel_${params.genome}.vcf.gz
    tabix -p vcf ${sample}_Nk_vep.vcf.gz
    """
}
