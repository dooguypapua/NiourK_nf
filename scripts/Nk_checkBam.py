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



#***** Check if BAM header correspund to reference FASTA
if len(sys.argv)!=5: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_checkBam refCompatibility] Missing arguments")
path_samtools = sys.argv[1]
if not os.path.isfile(path_samtools): exit("🅴 🆁 🆁 🅾 🆁\n[Nk_checkBam refCompatibility] Unable to find file `"+path_samtools+"`")
try: nb_threads = int(sys.argv[2])
except: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_checkBam refCompatibility] Invalid thread parameter `"+nb_threads+"`")
path_bam = sys.argv[3]
if not os.path.isfile(path_bam): exit("🅴 🆁 🆁 🅾 🆁\n[Nk_checkBam refCompatibility] Unable to find file `"+path_bam+"`")
path_ref_dict = sys.argv[4]
if not os.path.isfile(path_ref_dict): exit("🅴 🆁 🆁 🅾 🆁\n[Nk_checkBam refCompatibility] Unable to find file `"+path_ref_dict+"`")
dico_compare = {}
lst_diffsize = []
lst_absentInBam = []
lst_absentInDict = []
# Retrieve BAM header
process = subprocess.Popen([path_samtools+" view -@ "+str(nb_threads)+" -H "+path_bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
output, errors = process.communicate()
if process.returncode!=0: exit("🅴 🆁 🆁 🅾 🆁\n[Nk_checkBam refCompatibility] samtools view -H")
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
if errors!="": exit("🅴 🆁 🆁 🅾 🆁\n[Nk_checkBam refCompatibility]\n  Incompatible specified BAM & REF\n"+errors[:-1])