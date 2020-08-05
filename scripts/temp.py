import gffutils


excludeBiotype = ["pseudogene","C_region","J_segment","V_segment","V_segment_pseudogene","D_segment","J_segment_pseudogene","C_region_pseudogene","other"]
gffutilsDB = gffutils.FeatureDB("/home/dooguy/Dev_prog/Dev_Niourk/NiourK-nf/reference/GRCh37/GRCh37.gff.db", keep_order=True)
BED = open("/tmp/temp.bed",'w')
padding = 0

for gene in gffutilsDB.features_of_type("gene"):
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

