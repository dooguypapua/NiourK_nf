#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : Nk_checkTools.py                 |    
# description     : Check Tools Version.             |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20200521                         |
# version         : 0.1                              |
# python_version  : 3.8.2                            |
#==================================================================
# USAGE: Nk_checkTools.py toolName=toolPath ... toolName=toolPath |
# OUT  : Errors or JSON 
#==================================================================
import sys
import os
import json
import time
import subprocess


lstTools = ["gatk","deepvariant","vt","strelka","vep","vcfanno","samtools","bedtools","tvc","octopus","vcf-validator"]
dico_tools_version = { 'module':sys.argv[0],\
                       'date':time.strftime("%c %Z", time.localtime()) }
lstErrors = []
if len(sys.argv)==1:
    print("USAGE : python Nk_checkTools.py toolName=toolPath ... toolName=toolPath")
    lstErrors = ["Any input argument"]

for i in range(1,len(sys.argv),1):
    if not "=" in sys.argv[i]:
        lstErrors.append("Invalid input argument `"+sys.argv[i]+"`")
    else: 
        toolName = sys.argv[i].split("=")[0]
        toolPath = sys.argv[i].split("=")[1] # docker image version for deepvariant
        # Check if known tool name
        if not toolName in lstTools:
            lstErrors.append("Unknown tool name `"+toolName+"`")
        else:
            # File exist test already done in nextflow
            toolVersion = None
            # Check tool version
            if toolName=="gatk":
                try: toolVersion = subprocess.check_output(toolPath+" HaplotypeCaller --version",stderr=subprocess.STDOUT, shell=True).decode().split("Version: ")[1].split("\n")[0]
                except: lstErrors.append("Unable to check `"+toolName+"` version")
            elif toolName=="deepvariant":
                try: toolVersion = subprocess.check_output("docker images -a | grep \""+toolPath+"\" | sed -e s/\"\\s\\+\"/\"\\t\"/g | cut -f 2",stderr=subprocess.STDOUT, shell=True).decode().replace("\n","").replace("\n","")
                except: lstErrors.append("Unable to find `"+toolName+" `"+toolPath+"` docker version")
                toolPath = "docker"
                if toolVersion=="": lstErrors.append("Docker permission denied for `"+toolName+"`")
                if toolVersion.__contains__("Got permission denied"): lstErrors.append("Docker permission denied for `"+toolName+"`")
            elif toolName=="vt":
                try: toolVersion = subprocess.check_output(toolPath+" normalize -?",stderr=subprocess.STDOUT, shell=True).decode().split("normalize v")[1].split("\n")[0]  
                except: lstErrors.append("Unable to check `"+toolName+"` version")
            elif toolName=="strelka":
                try: toolVersion = subprocess.check_output(toolPath+" --version",stderr=subprocess.STDOUT, shell=True).decode().split("\n")[0]
                except: lstErrors.append("Unable to check `"+toolName+"` version")
            elif toolName=="vep":
                try: toolVersion = subprocess.check_output(toolPath+" | grep \"  ensembl\" | sed s/\" \"//g",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1]
                except: lstErrors.append("Unable to check `"+toolName+"` version")
            elif toolName=="vcfanno":
                try: toolVersion = subprocess.check_output(toolPath,stderr=subprocess.STDOUT, shell=True).decode().split("\n")[2].split(" ")[2]
                except: lstErrors.append("Unable to check `"+toolName+"` version")
            elif toolName=="vcf-validator":
                try: toolVersion = subprocess.check_output(toolPath+" -v",stderr=subprocess.STDOUT, shell=True).decode().split("\n")[0].split(" ")[2]
                except: lstErrors.append("Unable to check `"+toolName+"` version")
            else: # "samtools","bedtools","tvc","octopus"
                try: toolVersion = subprocess.check_output(toolPath+" --version",stderr=subprocess.STDOUT, shell=True).decode().split(" ")[1].split("\n")[0].replace("v","")
                except: lstErrors.append("Unable to check `"+toolName+"` version")
            # Add to dico
            if toolVersion: dico_tools_version[toolName] = {'path':toolPath, 'version':toolVersion}

# Write error if any
if len(lstErrors)>0: exit("\n".join(lstErrors))
else: print(json.dumps(dico_tools_version, indent=4))