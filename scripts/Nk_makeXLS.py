#=====================================================
# -*- coding: utf-8 -*-                               |
# title           : Nk_checkBed.py                    | 
# description     : Create a refseq region BED file   |
# author          : dooguypapua                       |
# copyright       : CHU Angers                        |
# date            : 20200809                          |
# version         : 0.1                               |
# python_version  : 3.8.2                             |
#==========================================================================================================
# USAGE: Nk_makeXLS.py pathVCF pathXLSX pathConf
# OUT  : 
#==========================================================================================================
import sys
import os
import xlsxwriter
import vcfpy
from Nk_functions import check_file_writable



#***** Arguments *****#
if len(sys.argv)<4: exit("[Nk_makeBED] Missing arguments")
if len(sys.argv)>4: exit("[Nk_makeBED] Too many arguments")
pathVCF = sys.argv[1]
if not os.path.isfile(pathVCF): exit("[Nk_makeXLS] Unable to find VCF file `"+pathBed+"`")
pathXLS = sys.argv[2]
if not check_file_writable(pathXLS):  exit("[Nk_makeXLS] Output XLS file not writable `"+pathXLS+"`")
pathConf = sys.argv[3]
if not os.path.isdir(pathConf) or os.path.basename(pathConf)!="conf": exit("[Nk_makeXLS] Invalid CONF folder `"+pathConf+"`")



#***** GET VEP variant consequence ontology *****#
dicoOntology = {}
ONTOLOGY = open(pathConf+"/vep_ontology.txt",'r')
lstLines = ONTOLOGY.read().split("\n")
ONTOLOGY.close()
for line in lstLines:
    if line!="" and line[0]!="#":
        splitLine = line.split("\t")
        dicoOntology[splitLine[0]] = splitLine[3]
if len(dicoOntology)==0: exit("[Nk_makeXLS] Unable to read VEP ontology file `"+pathConf+"/vep_ontology.txt`")
# Tools output reformat
dicoToolFormat = { \
                  "sift":{ "deleterious":"D+", "deleterious_low_confidence":"D", "tolerated":"B+", "tolerated_low_confidence":"B" } ,\
                  "polyphen":{ "probably_damaging":"D+", "possibly_damaging":"D", "benign":"B", "unknown":"?" } ,\
                  "clinvar":{ "uncertain_significance":"?", "not_provided":"?", "benign":"B+", "likely_benign":"B", "likely_pathogenic":"D", "pathogenic":"D+", "drug_response":"other", "histocompatibility":"other", "other":"other", "confers_sensitivity":"other", "risk_factor":"other", "association":"other", "protective":"other", "affects":"other" } \
                 }


#***** XLSXWRITER *****#
# Constructor
workbook = xlsxwriter.Workbook(pathXLS)
sheet = workbook.add_worksheet("NiourK")
cptRow = 0
cptCol = 0
# Format
cellFormat = workbook.add_format({'bold': 0,'border': 1,'align': 'center','valign': 'vcenter'})
headerFormat = workbook.add_format({'bold': 1,'border': 1,'align': 'center','valign': 'vcenter'})
hyperlinkFormat = workbook.add_format({'bold': 0,'underline':1,'border': 1,'align': 'center','valign': 'vcenter'})

dicoWidth = {}
# Header
lstHeader = ["CHR","POS","REF","ALT","VARSOME","NK FREQ","NK OCC","CSQ","GENE","TYPE","TRANSCRIPT","REGION","AF","DP","GT","QUAL","CALL","IMPACT","CLINVAR","GNOMAD","SIFT","PP2","CADD","SPLICEAI"]
for header in lstHeader:
    sheet.write(cptRow,cptCol,header,headerFormat) ; dicoWidth[cptCol] = len(header)*1.5
    cptCol+=1
cptRow+=1

# qual
# filter
# AF
# DP
# GT
# CALLAF
# CALLFILTER
# CALLNB
# CALLQUAL


#***** READ input VCF *****#
try:
    vcfReader = vcfpy.Reader.from_path(pathVCF)
    vcfHeader = vcfReader.header
except:
    exit("[Nk_makeXLS] Unable to read input VCF file `"+pathVCF+"`")
nkVersion = ""
lstCaller = []
dicoVepField = {}
# Retrieve header features
splitField = vcfHeader.get_info_field_info("CSQ").value.replace("\">","").split("Format: ")[1].split("|")
for i in range(len(splitField)): dicoVepField[splitField[i]] = i
for headerLine in vcfHeader.lines:
    if headerLine.key=="Nk_version": nkVersion = headerLine.value
    if headerLine.key=="Nk_calls": lstCaller.extend(headerLine.value.split("|"))
        
# Check well formated header
if nkVersion=="": exit("[Nk_makeXLS] Missing or empty `Nk_version` tag in input vcf header `"+pathVCF+"`")
if len(lstCaller)==0: exit("[Nk_makeXLS] Missing or empty `Nk_calls` tag in input vcf header `"+pathVCF+"`")
# Browse variant
for record in vcfReader:
    # Filter variants
    if record.FILTER[0]=="PASS":
        # Parse VEP annotation
        splitVEP = record.INFO["CSQ"][0].split(",")
        for annoVEP in splitVEP:
            cptCol = 0
            splitAnnoVep = annoVEP.split("|")
            for header in lstHeader:
                value = ""
                if header=="CHR": value = str(record.CHROM)
                elif header=="POS": value = int(record.POS)
                elif header=="REF": value = str(record.REF)
                elif header=="ALT": value = str(record.ALT[0].value)
                elif header=="QUAL": value = int(record.QUAL)
                elif header=="AF": value = float(record.calls[0].data.get('AF')[0])
                elif header=="DP": value = int(record.calls[0].data.get('DP'))
                elif header=="GT": value = str(record.calls[0].data.get('GT'))
                elif header=="CSQ":
                    for csq in splitAnnoVep[dicoVepField["Consequence"]].split("&"):
                        value+=str(dicoOntology[csq])+"&"
                    value = value[:-1].replace("&"," & ")
                elif header=="GENE": value = str(splitAnnoVep[dicoVepField["SYMBOL"]])
                elif header=="TYPE":
                    if splitAnnoVep[dicoVepField["BIOTYPE"]]=="protein_coding": value = "coding"
                    elif "pseudogene" in splitAnnoVep[dicoVepField["BIOTYPE"]]: value = "pseudo"
                    else: value = str(splitAnnoVep[dicoVepField["BIOTYPE"]])
                elif header=="TRANSCRIPT": value = str(splitAnnoVep[dicoVepField["Feature"]])
                elif header=="REGION":
                    if splitAnnoVep[dicoVepField["EXON"]]!="": value = "E"+splitAnnoVep[dicoVepField["EXON"]].split("/")[0]+" ("+splitAnnoVep[dicoVepField["EXON"]].split("/")[1]+")"
                    elif splitAnnoVep[dicoVepField["INTRON"]]!="": value = "I"+splitAnnoVep[dicoVepField["INTRON"]].split("/")[0]+" ("+splitAnnoVep[dicoVepField["INTRON"]].split("/")[1]+")"
                    else: value = ""
                elif header=="CALL": value = record.INFO["CALLNB"][0]
                elif header=="IMPACT": value = str(splitAnnoVep[dicoVepField["IMPACT"]]).title()
                elif header=="CLINVAR":
                    splitClinvar = set(splitAnnoVep[dicoVepField["CLIN_SIG"]].split("&"))
                    setValue = set()
                    for field in splitClinvar:
                        try: setValue.add(dicoToolFormat["clinvar"][field])
                        except: pass
                    value = " & ".join(list(setValue))
                elif header=="GNOMAD": value = float(splitAnnoVep[dicoVepField["gnomAD_AF"]])
                elif header=="SIFT":
                    try: value = str(dicoToolFormat["sift"][splitAnnoVep[dicoVepField["SIFT"]]])
                    except: value = str(splitAnnoVep[dicoVepField["SIFT"]])
                elif header=="PP2":
                    try: value = str(dicoToolFormat["polyphen"][splitAnnoVep[dicoVepField["PolyPhen"]]])
                    except: value = str(splitAnnoVep[dicoVepField["PolyPhen"]])
                elif header=="CADD": value = float(splitAnnoVep[dicoVepField["CADD_PHRED"]])
                elif header=="SPLICEAI":
                    if splitAnnoVep[dicoVepField["SpliceAI_pred_DS_AG"]]!="":
                        maxDeltaScore = 0
                        for spliceType in ["SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG","SpliceAI_pred_DS_DL"]:
                            if splitAnnoVep[dicoVepField[spliceType]]!="":
                                deltaScore = float(splitAnnoVep[dicoVepField[spliceType]])
                                if deltaScore>maxDeltaScore: maxDeltaScore = deltaScore
                        if maxDeltaScore>=0.8: value = "High"
                        elif maxDeltaScore>=0.5: value = "Moderate"
                        elif maxDeltaScore>=0.2: value = "Low"
                        else: value = "None"
                elif header=="VARSOME": value = "https://varsome.com/variant/hg38/"+str(record.CHROM)+"%3A"+str(record.POS)+"%3A"+str(record.REF)+"%3A"+str(record.ALT[0].value)
                if "http" in str(value):
                    sheet.write_url(cptRow,cptCol,value,hyperlinkFormat,string="link")
                else:
                    sheet.write(cptRow,cptCol,value,cellFormat)
                    strSize = sum(1 for c in str(value) if c.isupper())*1.2+sum(1 for c in str(value) if c.islower())+sum(1 for c in str(value) if c.isdigit())
                dicoWidth[cptCol] = max(strSize,dicoWidth[cptCol])
                cptCol+=1
            cptRow+=1


        # # Calling results
        # for i in range(len(lstCaller)):
        #     if record.INFO["CALLFILTER"][0].split("")[i]!=".":
        #         dicoVar["filter"][lstCaller[i]] = record.INFO["CALLFILTER"][0].split("")[i]
        #     else: 
        #         dicoVar["call"][lstCaller[i]] = record.INFO["CALLQUAL"][0].split("")[i]
        #     if record.INFO["CALLAF"][0].split("")[i]!=".":
        #         dicoVar["callaf"] = float(record.INFO["CALLAF"][0].split("")[i])
           

# Adjust Column width
for key in dicoWidth:
    sheet.set_column(key,key,dicoWidth[key]+1)

workbook.close()


# lstHeader = ["CSQ","GENE","TYPE","OMIM INH","TRANSCRIPT","LOCATION","CALL","IMPACT","GNOMAD"]
