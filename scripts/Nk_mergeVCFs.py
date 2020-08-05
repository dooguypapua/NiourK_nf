#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : Nk_mergeVCFs.py                  |    
# description     : Validate and Merge caller VCFs   |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20200525                         |
# version         : 0.2                              |
# python_version  : 3.8.2                            |
#==========================================================================================================
# USAGE: Nk_mergeVCFs.py path_gatk pathVT pathVcfValidator pathVCFraw pathFasta DPfilter sample VCFin ... VCFin
# OUT  : VCFmergeOut or Errors 
#==========================================================================================================
import sys
import os
import shutil
import subprocess
import vcfpy
import numpy
import collections


def sortVCF(pathVCFUnsorted,pathVCFSorted):
    cmd_sort_vcf = path_gatk+" SortVcf -I "+pathVCFUnsorted+" -O "+pathVCFSorted
    process = subprocess.Popen([cmd_sort_vcf], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] Sort merged VCF\n    "+err.decode('utf-8'))




def validateVCF(path_vcfvalidator,path_vcf):
    # Report folder
    path_dir_validate = path_vcf.replace(".vcf","_validateVCF")
    if os.path.isdir(path_dir_validate): shutil.rmtree(path_dir_validate)
    os.mkdir(path_dir_validate)
    # Create process
    cmd_validate = path_vcfvalidator+" -o "+path_dir_validate+" --require-evidence < "+path_vcf
    process = subprocess.Popen([cmd_validate], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    # retrieve statuts & summary file path
    statut = ""
    boolvalid = False
    reportPath = ""
    lst_errors = []
    for line in err.decode('utf-8').split("\n"):
        if line.__contains__("Summary report written to"): reportPath = line.split(": ")[1]
        elif line.__contains__("the input file is not valid"): statut = "not valid"
        elif line.__contains__("the input file is valid"): statut = "valid"
    if statut=="" or reportPath=="": lst_errors.append(["program error"])
    elif statut=="not valid":
        IN = open(reportPath,'r') ; lst_lines = IN.read().split("\n") ; IN.close()
        for line in lst_lines:
            if "Error:" in line: lst_errors.append(line.split("Error: ")[1])
    else: boolvalid = True
    shutil.rmtree(path_dir_validate)
    return (boolvalid,lst_errors)



#***** ARGUMENTS *****#
lst_vcf_arg = []
for i in range(1,len(sys.argv),1):
    if i==1: path_gatk = sys.argv[i]
    elif i==2: path_vt = sys.argv[i]
    elif i==3: path_vcfvalidator = sys.argv[i]
    elif i==4: pathVCFraw = sys.argv[i]
    elif i==5: pathFasta = sys.argv[i]
    elif i==6: DPfilter = sys.argv[i]
    elif i==7: sample = sys.argv[i]
    else: lst_vcf_arg.append(sys.argv[i])
# Check errors
errors = ""
# tools errors
if not os.path.isfile(path_gatk): errors+="    Unable to find GATK `"+path_gatk+"`\n"
process = subprocess.Popen([path_gatk], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: errors+="    Unable to launch GATK `"+path_gatk+"`\n"
elif not str(out).__contains__("GATK"): errors+="    Invalid GATK application `"+path_vt+"`\n"
if not os.path.isfile(path_vt): errors+="    Unable to find vt `"+path_vt+"`\n"
process = subprocess.Popen([path_vt+" -v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: errors+="    Unable to launch vt `"+path_vt+"`\n"
elif not str(err).__contains__("vt v"): errors+="    Invalid vt application `"+path_vt+"`\n"
if not os.path.isfile(path_vcfvalidator): errors+="    Unable to find vcf-validator `"+path_vcfvalidator+"`\n"
process = subprocess.Popen([path_vcfvalidator+" -v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: errors+="    Unable to launch vcf-validator `"+path_vcfvalidator+"`\n"
elif not str(out).__contains__("vcf_validator"): errors+="    Invalid vcf-validator application `"+path_vcfvalidator+"`\n"
# path errors
if not os.path.isdir(pathVCFraw):
    try: os.mkdir(pathVCFraw)
    except: errors+="    Unable to create output raw VCF folder `"+pathVCFraw+"`\n"
if not os.path.isfile(pathFasta): errors+="    Unable to find reference `"+pathFasta+"`\n"
elif not os.path.splitext(pathFasta)[1]==".fasta": errors+="    Invalid reference FASTA extension `"+pathFasta+"`\n"
try: dp_filter = int(DPfilter)
except:errors+="    Invalid DP filter value `"+DPfilter+"`\n"
lst_vcf_sample = []
for path_vcf in lst_vcf_arg: 
    if sample+"_" in os.path.basename(path_vcf):
        if not os.path.isfile(path_vcf): errors+="    Unable to find input VCF `"+path_vcf+"`\n"
        elif not os.path.splitext(path_vcf)[1]==".vcf": errors+="    Invalid input VCF extension `"+path_vcf+"`\n"
        else: lst_vcf_sample.append(path_vcf)
if len(path_vcf)==0: errors+="    Any correspunding VCF found for input sample `"+sample+"`\n"
if errors!="": exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] Arguments\n"+errors[:-1])


#***** VCF INPUT LOOP *****#
for path_vcf in lst_vcf_sample:

    path_temp_vcf = path_vcf.replace(".vcf","_formatencode.vcf")
    path_temp_temp_vcf = path_vcf.replace(".vcf","_formatencode_temp.vcf")
    path_filtered_vcf = path_vcf.replace(".vcf","_filtered.vcf")
    path_filteredSort_vcf = path_vcf.replace(".vcf","_filtered_sort.vcf")
    path_normalized_vcf = path_vcf.replace(".vcf","_normalize.vcf")

    #***** FORMAT & ENCODE *****#
    cmd_format_vcf = "sed -e 's/\\x0//g'" # 1 - replace bad tvc <0x00>
    cmd_format_vcf+= " -e 's/,assembly=[^,]*,/,/g'" # 2 - delete assembly info in TVC contig header
    cmd_format_vcf+= " -e '/##reference.*/d' -e 's/#CHROM/##reference="+pathFasta.replace("/","\\/")+"\\n#CHROM/g'" # 3 - add reference tag
    cmd_format_vcf+= " -e 's/AN=1/AN=2/g' "+path_vcf+" > "+path_temp_vcf # 4 - replace AN=1 to AN=2 for octopus
    cmd_format_vcf+= " ; old_encoding=$(file --brief --mime-encoding "+path_temp_vcf+")" # 5 - detect VCF encoding
    cmd_format_vcf+= " ; iconv -f $old_encoding -t utf-8 "+path_temp_vcf+" > "+path_temp_temp_vcf # 6 - convert encoding
    cmd_format_vcf+= " ; mv "+path_temp_temp_vcf+" "+path_temp_vcf
    process = subprocess.Popen([cmd_format_vcf], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] Format & encode\n    "+err.decode('utf-8'))

    #***** FILTER by DP *****#
    cmd_vcffilter = path_gatk+" VariantFiltration --output "+path_filtered_vcf+" --variant "+path_temp_vcf+" --filter-expression \"DP >= "+str(dp_filter)+"\" --filter-name \"DepthofQuality\""
    process = subprocess.Popen([cmd_vcffilter], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] DP VariantFiltration\n    "+err.decode('utf-8'))
    # reencode gatk date
    cmd_format_vcf+= "old_encoding=$(file --brief --mime-encoding "+path_filtered_vcf+")" # 5 - detect VCF encoding
    cmd_format_vcf+= " ; iconv -f $old_encoding -t utf-8 "+path_filtered_vcf+" > "+path_temp_temp_vcf # 6 - convert encoding
    cmd_format_vcf+= " ; mv "+path_temp_temp_vcf+" "+path_filtered_vcf

    #***** SORT VCF *****#
    sortVCF(path_filtered_vcf,path_filteredSort_vcf)

    #***** VALIDATE *****#
    boolvalid,lst_errors = validateVCF(path_vcfvalidator,path_filteredSort_vcf)
    if boolvalid==False: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] Validate VCF `"+os.path.basename(path_filteredSort_vcf)+"`\n    "+"\n    ".join(lst_errors))

    #***** COPY VCF to raw output folder *****#
    shutil.copy(path_filteredSort_vcf,pathVCFraw+"/"+os.path.basename(path_filteredSort_vcf))
    
    #***** DECOMPOSE & NORMALIZE *****#
    cmd_vt = path_vt+" decompose -s "+path_filteredSort_vcf+" | "+path_vt+" normalize -r "+pathFasta+" -o "+path_normalized_vcf+" -"
    process = subprocess.Popen([cmd_vt], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] Decompose & Normalize\n    "+err.decode('utf-8'))



#***** MERGE callers VCFs *****#
lst_caller_name = []
lst_contig_line = []
dico_filter_line = {}
dico_vcf = {}
pathMergeVCF = sample+"_Nk.vcf"
pathMergeUnsortedVCF = pathMergeVCF.replace(".vcf","_unsorted.vcf")
#***** INIT new vcf header *****#
new_header = vcfpy.Header(lines=None, samples=None)
new_header.add_line(vcfpy.HeaderLine("fileformat","VCFv4.2"))
#***** BROWSE caller vcf *****#
for path_vcf in lst_vcf_sample:
    caller_name = os.path.basename(path_vcf).split("_")[2].replace(".vcf","")
    lst_caller_name.append(caller_name)
    path_normalized_vcf = path_vcf.replace(".vcf","_normalize.vcf")
    vcf_tool_reader = vcfpy.Reader.from_path(path_normalized_vcf)
    vcf_header = vcf_tool_reader.header
    #***** READ HEADERS *****#
    # check header sample
    if new_header.samples==None: new_header.samples = vcf_header.samples
    # check header filters
    for filter_line in vcf_header.get_lines("FILTER"):
        if not filter_line.id in dico_filter_line: dico_filter_line[filter_line.id] = filter_line.description
        elif not filter_line.description in dico_filter_line[filter_line.id]: dico_filter_line[filter_line.id]+=", "+filter_line.description
    # check header contigs
    for contig_line in vcf_header.get_lines("contig"):
        if not contig_line in lst_contig_line: lst_contig_line.append(contig_line)
    #***** READ VARIANTS *****#
    for record in vcf_tool_reader:
        var_id = record.CHROM+"_"+str(record.POS)+"_"+record.REF+"_"+str(record.ALT[0].value)
        if not var_id in dico_vcf: dico_vcf[var_id] = { "features":{}, "ALT":record.ALT, "sample":new_header.samples.names[0]}
        dico_vcf[var_id]["features"][caller_name] = {}
        # Variant calling score (QUAL) field
        dico_vcf[var_id]["features"][caller_name]["QUAL"] = record.QUAL
        # Filter field
        dico_vcf[var_id]["features"][caller_name]["FILTER"] = record.FILTER
        # Genotype (GT) field
        dico_vcf[var_id]["features"][caller_name]["GT"] = record.calls[0].data.get('GT').replace("|","/")
        # Read Depth (DP) field
        dico_vcf[var_id]["features"][caller_name]["DP"] = record.calls[0].data.get('DP')
        # Allele Frequency (AF) field
        if caller_name in ["strelka","deepvariant"]: dico_vcf[var_id]["features"][caller_name]["AF"] = round((float(record.calls[0].data['AD'][1])/(float(record.calls[0].data['AD'][0])+float(record.calls[0].data['AD'][1]))),2) # for strelka and deepvariant AD for ref and alt is in FORMAT
        else: dico_vcf[var_id]["features"][caller_name]["AF"] = round(float(record.INFO['AF'][0]),2)
#***** CREATE new vcf header *****#
# Callers list
new_header.add_line(vcfpy.HeaderLine("Nk_calls","|".join(lst_caller_name)))
# Filters list
for filter_id in dico_filter_line: new_header.add_filter_line(vcfpy.OrderedDict([('ID', filter_id),('Description', dico_filter_line[filter_id])]))
new_header.add_filter_line(vcfpy.OrderedDict([('ID', "FILTER"),('Description', "All callers filtered")]))
# INFO
dico_line_info_CALLNB = collections.OrderedDict([("ID","CALLNB"),("Number","A"),("Type","Integer"),("Description","Number of PASS calls")])
new_header.add_info_line(dico_line_info_CALLNB)
dico_line_info_CALLAF = collections.OrderedDict([("ID","CALLAF"),("Number","A"),("Type","String"),("Description","Allele frequency per caller")])
new_header.add_info_line(dico_line_info_CALLAF)
dico_line_info_CALLFILTER = collections.OrderedDict([("ID","CALLFILTER"),("Number","A"),("Type","String"),("Description","Filters per caller")])
new_header.add_info_line(dico_line_info_CALLFILTER)
dico_line_info_CALLQUAL = collections.OrderedDict([("ID","CALLQUAL"),("Number","A"),("Type","String"),("Description","Variant quality per caller")])
new_header.add_info_line(dico_line_info_CALLQUAL)
# FORMAT
dico_line_info_GT = collections.OrderedDict([("ID","GT"),("Number","1"),("Type","String"),("Description","Genotype")])
new_header.add_format_line(dico_line_info_GT)
dico_line_info_DP = collections.OrderedDict([("ID","DP"),("Number","1"),("Type","Integer"),("Description","Read depth (median if multiple calls)")])
new_header.add_format_line(dico_line_info_DP)
dico_line_info_AF = collections.OrderedDict([("ID","AF"),("Number","A"),("Type","Float"),("Description","Allele Frequency (median if multiple calls)")])
new_header.add_format_line(dico_line_info_AF)
for contig_line in lst_contig_line: new_header.add_line(contig_line)
# REFERENCE
new_header.add_line(vcfpy.HeaderLine("reference",pathFasta))
# Write header
writer = vcfpy.Writer.from_path(pathMergeUnsortedVCF, new_header)
#***** CREATE new vcf variant lines *****#
for var_id in dico_vcf:
    # Main
    split_var_id = var_id.split("_")
    chrom = split_var_id[0]
    pos = int(split_var_id[1])
    ref = split_var_id[2]
    # Compute others fields
    set_gt = set()
    lst_dp = []
    lst_af = []
    lst_qual = []
    lst_filter = []
    nb_call = 0
    for caller_name in lst_caller_name:
        if caller_name in dico_vcf[var_id]["features"]:
            set_gt.add(dico_vcf[var_id]["features"][caller_name]["GT"])
            if dico_vcf[var_id]["features"][caller_name]["DP"]==None: lst_dp.append(0)
            else: lst_dp.append(dico_vcf[var_id]["features"][caller_name]["DP"])
            lst_af.append(dico_vcf[var_id]["features"][caller_name]["AF"])
            lst_qual.append(dico_vcf[var_id]["features"][caller_name]["QUAL"])
            lst_filter.append(" ".join(dico_vcf[var_id]["features"][caller_name]["FILTER"]))
            if "PASS" in dico_vcf[var_id]["features"][caller_name]["FILTER"]: nb_call+=1
        else:
            lst_af.append(".")
            lst_filter.append(".")
    # Qual
    field_qual = int(round(numpy.median(lst_qual),0))
    # Filter
    if "PASS" in lst_filter: field_filter = "PASS"
    else: field_filter = "FILTER"
    # Info
    dico_info = { "CALLNB":[nb_call], "CALLAF":["|".join(numpy.array(lst_af,dtype=str))], "CALLFILTER":["|".join(lst_filter).replace(" ","")], "CALLQUAL":["|".join(numpy.array(lst_qual,dtype=str))] }
    #***** FORMAT *****#
    lst_format_id = ['GT', 'DP', 'AF']
    # Merge GT field
    try: set_gt.remove('./.')
    except: pass
    if len(set_gt)==1: field_gt = set_gt.pop()
    else: field_gt = "./."
    # Merge DP field
    field_dp = int(round(numpy.median(lst_dp),0))
    # Merge AF field
    while lst_af.count(".")>0: lst_af.remove(".")
    field_af = float(round(numpy.median(lst_af),2))
    # Create call
    dico_calls = [vcfpy.Call(dico_vcf[var_id]["sample"], {'GT':field_gt, 'DP':field_dp, 'AF':[field_af]})] 
    #***** WRITE VARIANT *****#
    new_record = vcfpy.Record(chrom, pos, ".", ref, dico_vcf[var_id]["ALT"], field_qual, [field_filter], dico_info, lst_format_id, dico_calls)
    writer.write_record(new_record)
writer.close()



#***** POST-PROCESSING *****#
# Sort
sortVCF(pathMergeUnsortedVCF,pathMergeVCF)
# Validate
boolvalid,lst_errors = validateVCF(path_vcfvalidator,pathMergeVCF)
if boolvalid==False: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] Validate VCF `"+os.path.basename(pathMergeVCF)+"`\n    "+"\n    ".join(lst_errors))
# bgzip
cmd_bgzip = "bgzip -f "+pathMergeVCF
process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] bgzip VCF `"+os.path.basename(pathMergeVCF)+"`\n")
# tabix
cmd_tabix = "tabix -p vcf "+pathMergeVCF+".gz"
process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_mergeVCFs] tabix VCF `"+os.path.basename(pathMergeVCF)+".gz`\n")