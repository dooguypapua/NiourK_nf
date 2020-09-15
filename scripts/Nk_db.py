#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : Nk_db.py                         |    
# description     : NiourK-DB management             |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20200910                         |
# version         : 0.2                              |
# python_version  : 3.8.2                            |
#==================================================================
# USAGE: python Nk_db.py <subcommand> [options] [parameters]
#
#   Nk_db sub-commands include:
#
#     add-run         Add sequencing run from a parameter file
#     add-depth       Add sample depth from a BED file
#     add-vcf         Add variants from a VCF file
#     add-nksample    Add variants from a Nk_sample JSON file
#     --------------------------------------------------------
#     del-run         Delete sequencing run
#     del-sample      Delete sample depth and variants
#     --------------------------------------------------------
#     list-run        List sequencing runs
#     list-sample     List samples for a sequencing run
#     --------------------------------------------------------
#     get-variant     Get variant features
#     --------------------------------------------------------
#     stats           Print some database statistics
#     export          Export database to a specified file
#     --------------------------------------------------------
#     -h  --help      Print this help menu
#     -v  --version   Print tool version
#     -q  --quiet     Disable verbose output
#==================================================================
# BATCH:
#  ADD-RUN (nk_run & nk_sample)
#  bash add_batch.sh run <parameters_json folder> <threads>
#  ADD-DEPTH
#  bash add_batch.sh depth <GRCh38 all run BED folder> <threads>
#  ADD-NKSAMPLE
#  bash add_batch.sh var <nksample all run folder> <threads>
#==================================================================

import sys
import os
import shutil
import tempfile
from Nk_functions import *
from Nk_mongo import *



#***** INITIALIZATION *****#
# Main init dictionnary
dicoInit = {
            'white':"235;235;235", 'grey1':"200;200;200", 'grey2':"150;150;150",'blue1':"135;135;222",'blue2':"175;175;233", 'red':"255;85;85",'green':"113;180;120" ,\
            'colorBool':True, "quiet" : False, \
            'pathSrc':os.path.dirname(os.path.abspath("__file__")), 'pathDirTmp':tempfile.mkdtemp(), \
            'pathLiftover':os.path.dirname(os.path.abspath("__file__"))+"/liftOver", 'pathLiftChain':os.path.dirname(os.path.abspath("__file__"))+"/hg19ToHg38.over.chain.gz", \
            'mongoHost':"localhost", "mongoPort":"27018", 'nbChunk':25000, 'maxSevSelDelay':10 , 'truncateWidth':(30,50) \
           }
# MongoDB (sudo mongod --port 27018 --dbpath /media/dooguy/ultima_thule/niourkdb)
connectMongo(dicoInit)

# Arguments
NkDBargManager(sys.argv,dicoInit)

# Check NkDepth
# if dicoInit["collect_nk_depth"].count()!=3088286401: initNkDepth(dicoInit,dicoChrSize)



#***** DB CONSULTATION *****#
# List NkRun
if dicoInit["subCmd"]=="list-run": listRun(dicoInit)
# List NkSample
if dicoInit["subCmd"]=="list-sample": listSample(dicoInit)
# DB Summary
if dicoInit["subCmd"]=="stats": dbStats(dicoInit)
# Variant Summary
if dicoInit["subCmd"]=="get-variant": getVariant(dicoInit)



#***** ADD FUNCTIONS  *****#
# Add Run
if dicoInit["subCmd"]=="add-run": addRun(dicoInit)
# Add VCF
if dicoInit["subCmd"]=="add-vcf": addVcf(dicoInit)
# Add nk_sample
if dicoInit["subCmd"]=="add-nksample": addNkSample(dicoInit)
# Add depth
if dicoInit["subCmd"]=="add-depth": addDepth(dicoInit)



#***** POSTPROCESSING  *****#
# Clean temporary folder
shutil.rmtree(dicoInit["pathDirTmp"])
# Exit
exit("\n")