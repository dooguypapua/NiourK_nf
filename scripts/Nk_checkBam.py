#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : Nk_checkBam.py                   |    
# description     : Check BAM features               |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20200522                         |
# version         : 0.1                              |
# python_version  : 3.8.2                            |
#==================================================================
# USAGE: Nk_checkBam.py refCompatibility 
# OUT  :  
#==================================================================
import sys
import os
import re
import subprocess
import matplotlib
import matplotlib.pyplot as plt



#***** Check if BAM header correspund to reference FASTA
if sys.argv[1] == "refCompatibility":
    if len(sys.argv)!=6: exit("[Nk_checkBam refCompatibility] Missing arguments")
    path_samtools = sys.argv[2]
    if not os.path.isfile(path_samtools): exit("[Nk_checkBam refCompatibility] Unable to find file `"+path_samtools+"`")
    try: nb_threads = int(sys.argv[3])
    except: exit("[Nk_checkBam refCompatibility] Invalid thread parameter `"+nb_threads+"`")
    path_bam = sys.argv[4]
    if not os.path.isfile(path_bam): exit("[Nk_checkBam refCompatibility] Unable to find file `"+path_bam+"`")
    path_ref_dict = sys.argv[5]
    if not os.path.isfile(path_ref_dict): exit("[Nk_checkBam refCompatibility] Unable to find file `"+path_ref_dict+"`")
    dico_compare = {}
    lst_diffsize = []
    lst_absentInBam = []
    lst_absentInDict = []
    # Retrieve BAM header
    process = subprocess.Popen([path_samtools+" view -@ "+str(nb_threads)+" -H "+path_bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = process.communicate()
    if process.returncode!=0: exit("[Nk_checkBam refCompatibility] samtools view -H")
    for line in output.decode('utf-8').split("\n"):
        if line.__contains__("@SQ"):
            split_line = line.split("\t")
            chrom_name = split_line[1].replace("SN:","")
            chrom_size = int(split_line[2].replace("LN:",""))
            dico_compare[chrom_name] = chrom_size
    # Read reference dict
    DICT = open(path_ref_dict,'r')
    lst_lines = DICT.read().split("\n")
    DICT.close()
    lst_contig_dict = []
    for line in lst_lines:
        if line.__contains__("@SQ"):
            split_line = line.split("\t")
            chrom_name = split_line[1].replace("SN:","")
            chrom_size = int(split_line[2].replace("LN:",""))
            if not chrom_name in dico_compare: lst_absentInBam.append(chrom_name)
            elif dico_compare[chrom_name]!=chrom_size: lst_diffsize.append(chrom_name)
            lst_contig_dict.append(chrom_name)
    for contig in dico_compare:
        if not contig in lst_contig_dict: lst_absentInDict.append(chrom_name)
    # Check errors
    errors = ""
    if len(lst_absentInBam)>0: errors+="    "+str(len(lst_absentInBam)).rjust(3)+" contig(s) not found in BAM header\n"
    if len(lst_absentInDict)>0: errors+="    "+str(len(lst_absentInDict)).rjust(3)+" contig(s) not found in REF dict\n"
    if len(lst_diffsize)>0: errors+="    "+str(len(lst_diffsize)).rjust(3)+" contig(s) with ≠ length\n"
    if errors!="": exit("[Nk_checkBam refCompatibility]\n  Incompatible specified BAM & REF\n"+errors[:-1])



#***** PLOT BamIndexStats files *****#
elif sys.argv[1] == "plotBamidxstats":
    # Arguments
    if len(sys.argv)!=3: exit("[Nk_checkBam plotBamidxstats] Missing arguments")
    if not os.path.isfile(sys.argv[2]): exit("[Nk_checkBam plotBamidxstats] Unable to find file `"+sys.argv[2]+"`")
    path_bamidxstats = sys.argv[2]
    dico_count = {}
    #***** PARSE BamIndexStats merge file *****#
    IN = open(path_bamidxstats,'r')
    lst_results = IN.read().split("==>")
    IN.close()
    for results in lst_results:
        if results.__contains__(".bamidxstats <=="):
            filename = results.split(" ")[1].replace(".bamidxstats","")
            dico_count[filename] = {}
            nb_reads_all = 0
            lst_lines = results.split("\n")
            for line in lst_lines:
                if line.__contains__("NoCoordinateCount"): dico_count[filename]["unmapped"] = int(line.replace("NoCoordinateCount= ",""))
                search_nbreads = re.search("(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+",line)
                if search_nbreads:
                    dico_count[filename][search_nbreads.group(1)] = int(search_nbreads.group(3))
                    nb_reads_all+=int(search_nbreads.group(3))
            dico_count[filename]["all"] = nb_reads_all
    # Write error if any

    #***** MAKE PLOT *****#
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    for filename in sorted(dico_count.keys()):
        lst_rpkm = []
        for chrom in dico_count[filename]:
            if chrom!="all": lst_rpkm.append(float(dico_count[filename][chrom])/float(dico_count[filename]["all"]))
        ax.plot(lst_rpkm, label=filename)
    ax.set_xticks(range(len(dico_count[filename])-1))
    ax.set_xticklabels(list(dico_count[filename].keys())[:-1], rotation='vertical')
    ax.get_yaxis().set_visible(False)
    ax.xaxis.grid(True)
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0 * 0.2, box.y0*1.4, box.width*0.99, box.height*0.99])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
    plt.title("Reads number ratio per chomosome", fontsize=14, fontweight='bold')
    plt.savefig(path_bamidxstats.replace(".txt","Plot.png"), dpi=150)
