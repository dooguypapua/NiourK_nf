#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : Nk_checkBed.py                   |    
# description     : Create a refseq region BED file  |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20202906                         |
# version         : 0.2                              |
# python_version  : 3.8.2                            |
#==========================================================================================================
# USAGE: Nk_checkBed.py pathBed pathGenes pathGff pathOut padding
# OUT  : 
#==========================================================================================================
import sys
import os
import subprocess
import gffutils
from pybedtools import BedTool,cleanup



#***** Arguments *****#
if len(sys.argv)!=6: exit("[Nk_AnnotBed] Missing arguments")
pathBed = sys.argv[1]
if pathBed!="None" and not os.path.isfile(pathBed): exit("[Nk_AnnotBed] Unable to find bed file `"+pathBed+"`")
pathGenes = sys.argv[2]
if pathGenes!="None" and not os.path.isfile(pathGenes): exit("[Nk_AnnotBed] Unable to find genes list file `"+pathGenes+"`")
pathGff = sys.argv[3]
if not os.path.isfile(pathGff): exit("[Nk_AnnotBed] Unable to find gff file `"+pathGff+"`")
if not os.path.isfile(pathGff+".db"): exit("[Nk_AnnotBed] Unable to find gff.db file `"+pathGff+".db`")
pathBedDirOut = sys.argv[4]
if not os.path.isdir(pathBedDirOut): exit("[Nk_AnnotBed] Invalid bed output folder `"+pathBedDirOut+"`")
try: padding = int(sys.argv[5])
except: exit("[Nk_AnnotBed] Invalid padding parameter `"+padding+"`")
# Load GFF database
try: gffutilsDB = gffutils.FeatureDB(pathGff+".db", keep_order=True)
except: exit("[Nk_AnnotBed] Unable to load gff.db file `"+pathGff+".db`")



#***** RETRIEVE GENE DB OBJECT *****#
setGene = set()
# 1sT case: no bed & no genes => keep all genes
if pathBed=="None" and pathGenes=="None":
    for gene in gffutilsDB.features_of_type("gene"):
        setGene.add(gene)
# 2nd case: no bed & genes => search correspunding gene identifier
elif pathGenes!="None":
    # Read genes list file
    IN = open(pathGenes,'r')
    lstGenes = IN.read().split("\n")
    IN.close()
    if len(lstGenes)==0: exit("[Nk_AnnotBed] Any gene found in genes list file `"+pathGenes+"`")
    # Search genes id
    lstMissingGene = lstGenes
    for gene in gffutilsDB.features_of_type("gene"):
        if gene.attributes["Name"][0] in lstGenes:
            setGene.add(gene)
            lstMissingGene.remove(gene.attributes["Name"][0])
    # Error if any gene not found
    if len(lstMissingGene)>0: exit("[Nk_AnnotBed] Unable to find following gene(s): `"+",".join(lstMissingGene)+"`")
#3nd case: bed & no genes => search intersected gene identifier
else:
    bed = BedTool(pathBed)
    genes = BedTool(pathGff)
    # Search gff exons intersection
    for intersect_elem in genes+bed:
        if intersect_elem.fields[2]=="exon":
            exon = gffutilsDB[intersect_elem.attrs["ID"]]
            # retrieve correspunding transcript
            for rna in gffutilsDB.parents(exon, order_by='start'):
                for gene in gffutilsDB.parents(rna, featuretype='gene', order_by='start'):
                    setGene.add(gene)
     # delete created temp file
    cleanup(remove_all=True)



#***** CONSTRUCT ANNOTATED REGIONS BED FILE *****#
BED = open("/tmp/tempNiourk.bed",'w')
excludeBiotype = ["pseudogene","C_region","J_segment","V_segment","V_segment_pseudogene","D_segment","J_segment_pseudogene","C_region_pseudogene","other"]
for gene in setGene:
    # Mandatory field
    geneName = gene.attributes["Name"][0]
    hgnc = "."
    mirbase = None
    for dbxref in gene.attributes["Dbxref"]:
        if "HGNC" in dbxref: hgnc = dbxref.split(":")[-1]
        if "miRBase" in dbxref: mirbase = dbxref.split(":")[-1]
    biotype = gene.attributes["gene_biotype"][0]
    if not biotype in excludeBiotype:
        #***** BROWSE Transcript *****#
        for rna in gffutilsDB.children(gene, featuretype=("snRNA","RNase_P_RNA","tRNA","RNase_MRP_RNA","antisense_RNA","Y_RNA","telomerase_RNA","miRNA","vault_RNA","ncRNA","mRNA","snoRNA","lnc_RNA","rRNA","SRP_RNA"), order_by='start'):
            # Count exons
            if rna.strand=="+": exonCpt = 1
            else: exonCpt = len(list(gffutilsDB.children(rna, featuretype='exon', order_by='start')))
            #***** BROWSE Exons *****#
            for exon in gffutilsDB.children(rna, featuretype='exon', order_by='start'):
                # No transcript for tRNA  => use product (ex: tRNA-Ile)
                # No transcript for miRNA => use gene name (ex: MIR6859-1)
                # No transcript for rRNA (RNR1 and RNR2) => use gene name (ex: RNR1)
                if mirbase!=None: transcriptId = mirbase
                elif rna.attributes["gbkey"][0]=="tRNA": transcriptId = rna.attributes["product"][0]
                elif geneName in ["RNR1","RNR2"]: transcriptId = geneName
                else: transcriptId = rna.attributes["transcript_id"][0]
                # Split exon including CDS
                try:
                    lstcds = list(gffutilsDB.children(rna, featuretype='CDS', order_by='start'))
                    cdsStart = lstcds[0].start
                    cdsEnd = lstcds[-1].end
                    if exon.start<cdsStart:
                        if exon.end>cdsStart:
                            if cdsEnd<exon.end: # case of one exon including all the cds
                                if exon.strand=="+":
                                    BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(cdsStart-1+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR5pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                    BED.write(exon.chrom+"\t"+str(cdsEnd-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR3pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                else:
                                    BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(cdsStart-1+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR3pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                    BED.write(exon.chrom+"\t"+str(cdsEnd-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR5pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                BED.write(exon.chrom+"\t"+str(cdsStart-1-padding)+"\t"+str(cdsEnd+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                            else: 
                                if exon.strand=="+": BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(cdsStart-1+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR5pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                else: BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(cdsStart-1+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR3pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                BED.write(exon.chrom+"\t"+str(cdsStart-1-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                        elif exon.strand=="+": BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR5pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                        else: BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR3pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                    elif exon.end>cdsEnd:
                        if exon.start<cdsEnd:
                            if exon.strand=="+": BED.write(exon.chrom+"\t"+str(cdsEnd-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR3pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                            else: BED.write(exon.chrom+"\t"+str(cdsEnd-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR5pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                            BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(cdsEnd+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                        elif exon.strand=="+": BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR3pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                        else: BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"_UTR5pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                    else:
                        BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                except:                    
                    BED.write(exon.chrom+"\t"+str(exon.start-1-padding)+"\t"+str(exon.end+padding)+"\t"+geneName+"#"+hgnc+"#"+biotype+"#"+transcriptId+"#E"+str(exonCpt).zfill(2)+"pad"+str(padding)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                # Increment exon number
                if rna.strand=="+": exonCpt+=1
                else: exonCpt-=1
BED.close()



#***** SORT & COMPRESS & TABIX *****#
# Sort
cmd_sort_bed = "sort -k 1V,1 -k 2n,2 /tmp/tempNiourk.bed > "+pathBedDirOut+"/NiourK.bed"
process = subprocess.Popen([cmd_sort_bed], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: exit("[Nk_AnnotBed] Sort NiourK BED file\n    "+err.decode('utf-8'))
# Compress
cmd_bgzip_bed = "bgzip -c "+pathBedDirOut+"/NiourK.bed > "+pathBedDirOut+"/NiourK.bed.gz"
process = subprocess.Popen([cmd_bgzip_bed], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: exit("[Nk_AnnotBed] Bgzip NiourK BED file\n    "+err.decode('utf-8'))
# Tabix
cmd_tabix_bed = "tabix -p bed "+pathBedDirOut+"/NiourK.bed.gz"
process = subprocess.Popen([cmd_tabix_bed], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = process.communicate()
if process.returncode!=0: exit("[Nk_AnnotBed] Index NiourK BED file\n    "+err.decode('utf-8'))
# Clean
os.remove("/tmp/tempNiourk.bed")