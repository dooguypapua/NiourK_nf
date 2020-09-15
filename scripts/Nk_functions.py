#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : Nk_checkBed.py                   |    
# description     : Create a refseq region BED file  |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20200910                         |
# version         : 0.1                              |
# python_version  : 3.8.2                            |
#=====================================================
import os
import re
import json
import shutil
import math
from io import StringIO





#---------------------------------------------------------------#
#---------------------------------------------------------------#
#                       DISPLAY FUNCTIONS                       #
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#***** COLOR *****#
def printcolor(text,style,fg_rgb,bg_rgb,colorBool):
    if colorBool==True:
        if bg_rgb: print("\x1b["+style+";38;2;"+fg_rgb+";48;2;"+bg_rgb+"m"+text+"\x1b[0m",end="", flush=True)
        else: print("\x1b["+style+";38;2;"+fg_rgb+"m"+text+"\x1b[0m",end="", flush=True)
    else: print(text,end="", flush=True)

#***** DISPLAY ERRORS *****#
def display_errors(lst_error,dicoInit):
    if len(lst_error)>0:
        printcolor(" "+lst_error[0]+"\n","0",dicoInit['red'],None,dicoInit['colorBool'])
        for i in range(1,len(lst_error),1): printcolor("         "+lst_error[i]+"\n","0",dicoInit['red'],None,dicoInit['colorBool'])
        try: shutil.rmtree(dicoInit["pathDirTmp"])
        except: pass
        exit("\n")

#***** CONVERT BYTE FORMAT *****#
def convertByteSize(size_bytes):
   if size_bytes == 0:
       return "0B"
   size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return "%s %s" % (s, size_name[i])





#---------------------------------------------------------------#
#---------------------------------------------------------------#
#                    HELP & USAGE FUNCTIONS                     #
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#***** HEADER *****#
def NkDBheaderUsage(dicoInit,boolUsage):
    printcolor("\nＮｉｏｕｒＫ • ＤＢ\n","1",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("  Python toolset\n\n","3",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("Version:     v0.2\n","0",dicoInit['grey2'],None,dicoInit['colorBool'])
    printcolor("Copyright:   CHU Angers\n","0",dicoInit['grey2'],None,dicoInit['colorBool'])
    printcolor("Developer:   David Goudenege\n","0",dicoInit['grey2'],None,dicoInit['colorBool'])
    printcolor("Code:        https://github.com/\n\n","0",dicoInit['grey2'],None,dicoInit['colorBool'])
    if boolUsage: printcolor("USAGE: ","1",dicoInit['blue1'],None,dicoInit['colorBool'])

#***** FOOTER *****#
def NkDBfooterUsage(dicoInit,error):
    printcolor("    --------------------------------------------------------\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    -h  --help      Print this help menu\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    -v  --version   Print tool version\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    -q  --quiet     Disable verbose output\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("\n"+error+"\n\n","0",dicoInit['red'],None,dicoInit['colorBool'])
    exit()

#***** MAIN USAGE *****#
def NkDBmainUsage(dicoInit,error):
    NkDBheaderUsage(dicoInit,True)
    printcolor("python Nk_db.py <subcommand> [options] [parameters]\n\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("  Nk_db sub-commands include:\n\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    add-run         Add sequencing run from a parameter file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    add-depth       Add sample depth from a BED file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    add-vcf         Add variants from a VCF file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    add-nksample    Add variants from a Nk_sample JSON file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    --------------------------------------------------------\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    del-run         Delete sequencing run\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    del-sample      Delete sample depth and variants\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    --------------------------------------------------------\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    list-run        List sequencing runs\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    list-sample     List samples for a sequencing run\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    --------------------------------------------------------\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    get-variant     Get variant features\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    --------------------------------------------------------\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    stats           Print some database statistics\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    export          Export database to a specified file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    NkDBfooterUsage(dicoInit,error)

#***** SUBCOMMAND USAGE *****#
# Add Run usage
def NkDBaddRunUsage(dicoInit,error):
    NkDBheaderUsage(dicoInit,True)
    printcolor("python Nk_db.py add-run --input <file>\n\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("    -i  --input     Input sequencing parameter file [required]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("                    (IonTorrent JSON or a Illumina XML)\n","3",dicoInit['grey1'],None,dicoInit['colorBool'])
    NkDBfooterUsage(dicoInit,error)
# Add Depth usage
def NkDBaddDepthNkSampleUsage(dicoInit,error):
    NkDBheaderUsage(dicoInit,True)
    printcolor("python Nk_db.py "+dicoInit["subCmd"]+" --input <file> --run <name/id> --sample <name/barcode>\n\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    if dicoInit["subCmd"]=="add-nksample": printcolor("    -i  --input     Input Nk_sample JSON file [required]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    else: printcolor("    -i  --input     Input sample depth BED(.gz) file [required]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    -r  --run       Sequencing name or id [required]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    -s  --sample    Sample name or barcode [required]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    NkDBfooterUsage(dicoInit,error)
# List Run & sample usage
def NkDBlistRunSampleUsage(dicoInit,error):
    NkDBheaderUsage(dicoInit,True)
    printcolor("python Nk_db.py "+dicoInit["subCmd"]+" --csv <file> [--pretty]"+"\n\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("    -o  --output    Output CSV file [required]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    -p  --pretty    Display summary [optionnal]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    NkDBfooterUsage(dicoInit,error)
# List Run & sample usage
def NkDBSimpleUsage(dicoInit,error):
    NkDBheaderUsage(dicoInit,True)
    printcolor("python Nk_db.py "+dicoInit["subCmd"]+"\n\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    NkDBfooterUsage(dicoInit,error)
# Get variant features
def NkDBaddGetVarUsage(dicoInit,error):
    NkDBheaderUsage(dicoInit,True)
    printcolor("python Nk_db.py get-variant --input <variant> [--mindepth <int>]\n\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("    -i  --input     Input variant [required]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("                    (must be formatted as follows: `chr5_145000789_A_G`)\n","3",dicoInit['grey1'],None,dicoInit['colorBool'])
    printcolor("    -m  --mindepth  Min depth to consider a sample covering a position [optionnal] [default:20]\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    NkDBfooterUsage(dicoInit,error)





#---------------------------------------------------------------#
#---------------------------------------------------------------#
#                  NIOURK-DB ARGUMENT MANAGER                   #
#---------------------------------------------------------------#
#---------------------------------------------------------------#
def NkDBargManager(lstArgv,dicoInit):

    #***** Help or no arguments => Main usage *****#
    if len(lstArgv)<2 or lstArgv[1] in ["--help","-h"]:
        NkDBmainUsage(dicoInit,"")

    #***** Subcommand *****#
    else:
        dicoInit["subCmd"] = lstArgv[1]
        
        #***** ADD run *****#
        if dicoInit["subCmd"]=="add-run":
            if len(lstArgv)<3 or lstArgv[2] in ["--help","-h"]: NkDBaddRunUsage(dicoInit,"")
            elif not lstArgv[2] in ["--input","-i"]: NkDBaddRunUsage(dicoInit,"Invalid argument `"+lstArgv[2]+"`")
            elif len(lstArgv)<4: NkDBaddRunUsage(dicoInit,"Missing value for `--input -i`")
            elif not os.path.isfile(lstArgv[3]): NkDBaddRunUsage(dicoInit,"Input file not found `"+lstArgv[3]+"`")
            else: dicoInit["pathInput"] = lstArgv[3]
        
        #***** ADD depth *****#
        elif dicoInit["subCmd"] in ["add-depth","add-nksample","add-vcf"]:
            if len(lstArgv)<3 or lstArgv[2] in ["--help","-h"]: NkDBaddDepthNkSampleUsage(dicoInit,"")
            if len(set(["--input","-i"]) & set(lstArgv))==0: NkDBlistRunSampleUsage(dicoInit,"Missing argument `--input -i`")
            if len(set(["--run","-r"]) & set(lstArgv))==0: NkDBlistRunSampleUsage(dicoInit,"Missing argument `--run -r`")
            if len(set(["--sample","-s"]) & set(lstArgv))==0: NkDBlistRunSampleUsage(dicoInit,"Missing argument `--sample -s`")
            dicoInit["pathInput"] = ""
            dicoInit["runId"] = ""
            dicoInit["runName"] = ""
            dicoInit["sampleBc"] = ""
            dicoInit["sampleName"] = ""
            for i in range(2,len(lstArgv),2):
                if len(lstArgv)<=i+1: NkDBaddDepthNkSampleUsage(dicoInit,"Missing value for `"+lstArgv[i]+"`")
                if lstArgv[i] in ["--input","-i"]:
                    if not os.path.isfile(lstArgv[i+1]): NkDBaddDepthNkSampleUsage(dicoInit,"Input file not found `"+lstArgv[i+1]+"`")
                    else: dicoInit["pathInput"] = lstArgv[i+1]
                elif lstArgv[i] in ["--run","-r"]:
                    # Search run entry
                    findRun = dicoInit["collect_nk_run"].find_one({"_id": lstArgv[i+1]})
                    if findRun:
                        dicoInit["runId"] = lstArgv[i+1]
                        dicoInit["runName"] = findRun["name"]
                    else:
                        findRun = dicoInit["collect_nk_run"].find_one({"name": lstArgv[i+1]})
                        if findRun:
                            dicoInit["runName"] = lstArgv[i+1]
                            dicoInit["runId"] = findRun["_id"]
                    try:
                        findRun = dicoInit["collect_nk_run"].find_one({"name": "_".join(lstArgv[i+1].split("_")[:-1]), "num":int(lstArgv[i+1].split("_")[-1]) })
                        if findRun:
                            dicoInit["runName"] = "_".join(lstArgv[i+1].split("_")[:-1])
                            dicoInit["runId"] = findRun["_id"]
                    except: pass
                elif lstArgv[i] in ["--sample","-s"]:
                    sample = lstArgv[i+1]
            # Check
            if dicoInit["pathInput"]=="": NkDBaddDepthNkSampleUsage(dicoInit,"Missing value for `--input -i`")
            if dicoInit["runId"]==dicoInit["runName"]=="": NkDBaddDepthNkSampleUsage(dicoInit,"Sequencing name or id not found in database")
            # Search sample entry
            findSample = dicoInit["collect_nk_sample"].find_one({"runid":dicoInit["runId"], "bc": lstArgv[i+1]})
            if findSample:
                dicoInit["sampleBc"] = lstArgv[i+1]
                dicoInit["sampleName"] = findSample["name"]
            else:
                findSample = dicoInit["collect_nk_sample"].find_one({"runid":dicoInit["runId"],"name": lstArgv[i+1]})
                if findSample:
                    dicoInit["sampleName"] = lstArgv[i+1]
                    dicoInit["sampleBc"] = findSample["bc"]
            if dicoInit["sampleBc"]==dicoInit["sampleName"]=="": NkDBaddDepthNkSampleUsage(dicoInit,"Sample name or barcode not found in database")
        
        #***** DEL run or sample *****#
        elif dicoInit["subCmd"]=="del-run": exit("\nDEL RUN\n\n")
        elif dicoInit["subCmd"]=="del-sample": exit("\nDEL SAMPLE\n\n")

        #***** LIST runs or samples *****#
        elif dicoInit["subCmd"] in ["list-run","list-sample"]:
            if len(lstArgv)==3 and lstArgv[2] in ["--help","-h"]: NkDBlistRunSampleUsage(dicoInit,"")
            if len(set(["--output","-o"]) & set(lstArgv))==0: NkDBlistRunSampleUsage(dicoInit,"Missing argument `--output -o`")
            if len(set(["--pretty","-p"]) & set(lstArgv))==1: dicoInit["listPretty"] = True
            else: dicoInit["listPretty"] = False
            dicoInit["pathOutput"] = ""
            i = 2
            while i < len(lstArgv):
                if len(lstArgv)<=i+1 and not lstArgv[i] in ["-p","--pretty"]: NkDBlistRunSampleUsage(dicoInit,"Missing value for `"+lstArgv[i]+"`")
                if lstArgv[i] in ["--output","-o"]: dicoInit["pathOutput"] = lstArgv[i+1] ; i+=2
                elif lstArgv[i] in ["-p","--pretty"]: i+=1
                else: NkDBlistRunSampleUsage(dicoInit,"Unknwon optionnal argument for `list-sample`")
            # Check
            if dicoInit["pathOutput"] in ["","-p","--pretty"]: NkDBlistRunSampleUsage(dicoInit,"Missing value for `--output -o`")

        #***** GET VARIANT FEATURES *****#
        elif dicoInit["subCmd"]=="get-variant":
            if len(lstArgv)<4 or lstArgv[2] in ["--help","-h"]: NkDBaddGetVarUsage(dicoInit,"")
            if len(set(["--input","-i"]) & set(lstArgv))==0: NkDBaddGetVarUsage(dicoInit,"Missing argument `--input -i`")
            dicoInit["varInput"] = ""
            dicoInit["mindepth"] = 20
            for i in range(2,len(lstArgv),2):
                if len(lstArgv)<=i+1: NkDBaddGetVarUsage(dicoInit,"Missing value for `"+lstArgv[i]+"`")
                if lstArgv[i] in ["--input","-i"]: dicoInit["varInput"] = lstArgv[i+1]
                elif lstArgv[i] in ["--mindepth","-m"]:
                    try: dicoInit["mindepth"] = int(lstArgv[i+1])
                    except: NkDBaddGetVarUsage(dicoInit,"Bad integer value for `"+lstArgv[i]+"`")
            if dicoInit["varInput"]=="": addGetVarUsage(dicoInit,"Missing value for `--input -i`")
            elif not re.search("chr[123456789YXM][0-9]?_[0-9]+_[ATGC]+_[ATGC]+",lstArgv[3]): NkDBaddGetVarUsage(dicoInit,"Invalid variant format `"+lstArgv[3]+"`")

        #***** STATS *****#
        elif dicoInit["subCmd"]=="stats":
            if len(lstArgv)==3 and lstArgv[2] in ["--help","-h"]: NkDBSimpleUsage(dicoInit,"")
            if len(lstArgv)>=3: NkDBSimpleUsage(dicoInit,"Any arguments required for `stats`")

        #***** DISPLAY version or error *****#
        elif dicoInit["subCmd"]=="--version" or dicoInit["subCmd"]=="-v":
            exit("Nk_db.py v2.0")
        else:
            NkDBmainUsage(dicoInit,"Invalid subcommand `"+dicoInit["subCmd"]+"`")
        # Quiet mode
        if "--quiet" in lstArgv or "-q" in lstArgv: dicoInit["quiet"] = True
    if not dicoInit["quiet"]: NkDBheaderUsage(dicoInit,False)






#---------------------------------------------------------------#
#---------------------------------------------------------------#
#                         INPUT/OUTPUT                          #
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#***** LOAD IonTorrent JSON  *****#
def loadTorrentJson(pathJson):
    JSON=open(pathJson)
    dataJson = json.load(JSON)
    JSON.close()
    try: HashExp = json.loads(dataJson['exp_json'])
    except: HashExp = dataJson['exp_json']
    try : HashBarcode =  json.load(StringIO(dataJson["experimentAnalysisSettings"]["barcodedSamples"]))
    except: HashBarcode = dataJson["experimentAnalysisSettings"]["barcodedSamples"]
    if type(HashBarcode)==str: json.load(json.load(StringIO(dataJson["experimentAnalysisSettings"]["barcodedSamples"])))
    return dataJson,HashExp,HashBarcode

#***** CHECK if output path could be write *****#
def check_file_writable(fnm):
    if os.path.exists(fnm):
        # path exists
        if os.path.isfile(fnm): # is it a file or a dir?
            # also works when file is a link and the target is writable
            return os.access(fnm, os.W_OK)
        else:
            return False # path is a dir, so cannot write as a file
    # target does not exist, check perms on parent dir
    pdir = os.path.dirname(fnm)
    if not pdir: pdir = '.'
    # target is creatable if parent dir is writable
    return os.access(pdir, os.W_OK)


