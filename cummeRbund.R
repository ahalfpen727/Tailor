###################################################
### Library+Variable --> 1) Load 2) Initiation
###################################################
options(width=65)
source("http://bioconductor.org/biocLite.R")
#biocLite("cummeRbund")
library(cummeRbund)
#source('openGraphSaveGraph.R');

# causes a traceback to appear if there is an error, and the traceback has the line number prefixed by {#}
options(error=traceback)

#options(echo=TRUE)
#args = commandArgs(trailingOnly = T)

# Gets arguments that were passed in via command line
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

###################################################
### code chunk number 1: args
###################################################

# These arguments are passed in on the command line via launchR.eachDiffDir.sh
# args:
# diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"  under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\" FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\"


######################################################
### code chunk number 1.5: Read Environment Variables
######################################################

print(args)

geneListString<-geneList
geneListString

# Cannonical Comman
cuff<-readCufflinks(dir=inDir,gtfFile=refgtf,genome=genomeR,out=outDir,rebuild=T)

# the following command is for rerunning cummerbund after gene level failure occurred
#cuff<-readCufflinks(dir=inDir,out=outDir, rebuild=F)
cuff

pdf(file.path(Rplots))

##############################################################################################
# code chunk for debugging purposes
# Must have all cuffdiff output (from *-over-* directory)
# as well as ref genome fasta and annotation gtf file
# in the current working directory
#############################################################################################

# geneListString<-c("CXCL12","TGFB","MYC","RAS","CXCR4","IL8","IL6")
# cuff<-readCufflinks(dir=".", gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
# refgtf="cuffcmp.combined.gtf"
# genome_path="genome.fa"
# over="LUTS"
# under="CTRL"

###########################################################################
### Chunk_2: Significantly differentially expressed features: method 1
############################################################################

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
mySigGeneTable<-getSigTable(cuff,alpha=0.05,level='genes')
head(mySigGeneTable,20)
length(mySigGeneTable)
sigGenes<-getGenes(cuff, mySigGenes)
sigGenes
head(sigGenes)length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)
sig_genes_annot<-sigGenes@annotation
head(sig_genes_exp.diff)
head(sig_genes_annot)
sig_genes_exp.diff["gene_id"]<-sig_genes_annot["gene_short_name"]
head(sig_genes_exp.diff)


mySigCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
head(mySigCDS)
length(mySigCDS)
mySigCDSTable<-getSigTable(cuff,alpha=0.05,level='CDS')
head(mySigCDSTable,20)
length(mySigCDSTable)
sigCDS<-getGenes(cuff, mySigCDS)
sigCDS
length(sigCDS)

sig_cds_exp.diff<-diffData(sigCDS)
sig_cds_annot<-sigCDS@annotation
head(sig_cds_exp.diff)
head(sig_cds_annot)
sig_cds_exp.diff["gene_id"]<-sig_cds_annot["gene_short_name"]
head(sig_cds_exp.diff)


mySigIsos<-getSig(cuff,x=over,y=under,alpha=.05,level='isoforms')
head(mySigIsos)
length(mySigIsos)
mySigIsoTable<-getSigTable(cuff,alpha=0.05,level='isoforms')
head(mySigIsoTable,20)
length(mySigIsoTable)
sigIsos<-getGenes(cuff, mySigIsos)
sigIsos
length(sigIsos)

sig_isos_exp.diff<-diffData(sigIsos)
sig_isos_annot<-sigIsos@annotation
head(sig_isos_exp.diff)
head(sig_isos_annot)
sig_isos_exp.diff["gene_id"]<-sig_isos_annot["gene_short_name"]
head(sig_isos_exp.diff)
mySigTSS<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
head(mySigTSS)
length(mySigTSS)
mySigTSSTable<-getSigTable(cuff,alpha=0.05,level='TSS')
head(mySigTSSTable,20)
length(mySigTSSTable)
sigTSS<-getGenes(cuff, mySigTSS)
sigTSS
length(sigTSS)

sig_tss_exp.diff<-diffData(sigTSS)
sig_tss_annot<-sigTSS@annotation
head(sig_tss_exp.diff)
head(sig_tss_annot)
sig_tss_exp.diff["gene_id"]<-sig_tss_annot["gene_short_name"]
head(sig_tss_exp.diff)
dim(sig_tss_exp.diff)

###################################################
### Chunk_3: Differentially Loaded Features
###################################################

mySigrelCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='relCDS')
head(mySigrelCDS)
length(mySigrelCDS)
sigrelCDS<-getGenes(cuff, mySigrelCDS)
sigrelCDS
length(sigrelCDS)
relcds.diff<-distValues(relCDS(cuff))
relcds.sigdiff<-subset(relcds.diff, significant=="yes")
head(relcds.sigdiff)
dim(relcds.sigdiff)

mySigSplices<-getSig(cuff,x=over,y=under,alpha=.05,level='splicing')
head(mySigSplices)
length(mySigSplices)
sigSplices<-getGenes(cuff, mySigSplices)
sigSplices
length(sigSplices)
splice.diff<-distValues(splicing(cuff))
splice.sigdiff<-subset(splice.diff, significant=="yes")
head(splice.sigdiff)
dim(splice.sigdiff)

mySigPromoters<-getSig(cuff,x=over,y=under,alpha=.05,level='promoters')
head(mySigPromoters)
length(mySigPromoters)
sigPromoters<-getGenes(cuff, mySigPromoters)
sigPromoters
length(sigPromoters)
promoter.diff<-distValues(promoters(cuff))
promoter.sigdiff<-subset(promoter.diff, significant=="yes")
head(promoter.sigdiff)
dim(promoter.sigdiff)

###################################################
### Chunk_4: Sig Matrix_plots
###################################################

mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
print(mySigMat)
sigIsoformIds1<-sigMatrix(cuff, level='isoforms',alpha=0.05)
print(sigIsoformIds1)
sigIsoformIds2<-sigMatrix(cuff, level='TSS',alpha=0.05)
print(sigIsoformIds2)
sigIsoformIds3<-sigMatrix(cuff, level='CDS', alpha=0.05)
print(sigIsoformIds3)

###############################################################################
# Chunk_5: Features, Counts, and all inclusive tables
# (CDS-p_ids, TSS-TSS ids, isoforms- XLOC and TCONS ids, and Genes XLOC ids)
################################################################################

runInfo(cuff)
replicates(cuff)
conditions(cuff)

gene.counts<-count(genes(cuff))
head(gene.counts)
gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)

iso.counts<-count(isoforms(cuff))
head(iso.counts)
iso.featurenames<-featureNames(isoforms(cuff))
head(iso.featurenames)

cds.counts<-count(CDS(cuff))
head(cds.counts)
cds.featurenames<-featureNames(CDS(cuff))
head(cds.featurenames)

tss.counts<-count(TSS(cuff))
head(tss.counts)
tss.featurenames<-featureNames(TSS(cuff))
head(tss.featurenames)

###################################################################
### Chunk_7: get gene symbols for CDS TSS Isoforms and Genes
#####################################################################

gene.xloc.matrix<-featureNames(genes(cuff))
head(gene.xloc.matrix)
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)

isos.tcons.matrix<-featureNames(isoforms(cuff))
head(isos.tcons.matrix)
iso.list<-getGenes(cuff,geneId = isos.tcons.matrix)
iso.list
iso_annotation_data<-featureNames(iso.list)
head(iso_annotation_data)

cds.p_id.matrix<-featureNames(CDS(cuff))
head(cds.p_id.matrix)
cds.list<-getGenes(cuff,geneId = cds.p_id.matrix)
cds.list
cds_annotation_data<-featureNames(cds.list)
head(cds_annotation_data)

tss.tssgroup.matrix<-featureNames(TSS(cuff))
head(tss.tssgroup.matrix)
tss.list<-getGenes(cuff,geneId = tss.matrix)
tss.list
tss_annotation_data<-featureNames(tss.list)
head(tss_annotation_data)

####################################################################
### code chunk number 8: Write repFPKMMatrix and DiffTable to file
####################################################################

g.rep.matrix<-repFpkmMatrix(genes(cuff))
gene.xloc.matrix<-featureNames(genes(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list

gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)

gene_rep_matrix<-cbind2(gene_annotation_data, g.rep.matrix)
gene.rep.matrix<-gene_rep_matrix[,-1]

fpkm.file = file.path(FPKMmatrix)
write.table(gene.rep.matrix,file = fpkm.file, sep = "  ", row.names = F, col.names = T,quote = F)

# This table are huge and chew up a lot of memory
isodiff<-diffTable(isoforms(cuff))
diff.table = file.path(DiffTable)
write.table(isodiff, file = diff.table, sep = "  ", row.names = F , col.names = T,quote = F)


##################################################
### code chunk number 10: global_dispersion
###################################################

disp_cuff<-dispersionPlot(cuff) + ggtitle("Overloading of All Cuffdiff Output")
disp_cuff
disp_gene<-dispersionPlot(genes(cuff)) + ggtitle("CUFFDIFF Genes Overloading")
disp_gene

###################################################
### code chunk number 11: SCV_visualization
###################################################

genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
genes.scv_pool<-fpkmSCVPlot(genes(cuff), showPool=T)
genes.scv_pool

isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
isoforms.scv
isoforms.scv_pool<-fpkmSCVPlot(isoforms(cuff), showPool=T)
isoforms.scv_pool

cvtss<-fpkmSCVPlot(TSS(cuff))
cvtss
cvtss_pool<-fpkmSCVPlot(TSS(cuff), showPool=T)
cvtss_pool
cvcds1<-fpkmSCVPlot(CDS(cuff))
cvcds1
cvcds_pool<-fpkmSCVPlot(TSS(cuff), showPool=T)
cvcds_pool

###################################################
### code chunk number 12: global Density_plots_1
###################################################
dens2<-csDensity(genes(cuff),replicates=F)
print(dens2)
dens3<-csDensity(isoforms(cuff), replicates=F)
print(dens3)
dens4<-csDensity(CDS(cuff), replicates=F)
print(dens4)
dens6<-csDensity(TSS(cuff), replicates=F)
print(dens6)

###################################################
### code chunk number 12: Per Sample Density_plots_1
###################################################
dens2<-csDensity(genes(cuff),replicates=T)
print(dens2)
dens3<-csDensity(isoforms(cuff), replicates=T)
print(dens3)
dens4<-csDensity(CDS(cuff), replicates=T)
print(dens4)
dens6<-csDensity(TSS(cuff), replicates=T)
print(dens6)

###################################################
### code chunk number 13: global_plots_2
###################################################
box_g<-csBoxplot(genes(cuff),replicates=T, logMode=T) + ggtitle("genes")
box_g
box_i<-csBoxplot(isoforms(cuff),replicates=T, logMode=T) + ggtitle("isoforms")
box_i
box_c<-csBoxplot(CDS(cuff),replicates=T, logMode=T) + ggtitle("CDS")
box_c
box_t<-csBoxplot(TSS(cuff),replicates=T, logMode=T) + ggtitle("TSS")
box_t
###################################################
### code chunk number 14: global_plots_scatter_1
###################################################

gene_scatter<-csScatter(genes(cuff),over,under,smooth=T) + ggtitle("genes")
gene_scatter
iso_scatter<-csScatter(isoforms(cuff), over, under,smooth=T, logMode=T) + ggtitle("isoforms")
iso_scatter
tss_scatter<-csScatter(TSS(cuff), over, under,smooth=T) + ggtitle("TSS")
tss_scatter
cds_scatter<-csScatter(CDS(cuff), over, under,smooth=T) + ggtitle("CDS")
cds_scatter
scgene<-csScatter(sigGenes,over,under,smooth=T) + ggtitle("sig Genes")
scgene
sciso<-csScatter(sigIsos,over, under,smooth=T) + ggtitle("Sig Isoforms")
sciso
sctss<-csScatter(sigTSS,over, under,smooth=F) + ggtitle("Sig TSS Groups")
sctss
sctss<-csScatter(sigCDS,over, under,smooth=F) + ggtitle("Sig CDS")
sctss
##################################################
### code chunk number 15: global_plots_dendro
###################################################

csDendro(genes(cuff), replicates=T)
csDendro(isoforms(cuff),replicates=T,pseudocount = 1, logMode=T)
csDendro(TSS(cuff),replicates=T)
csDendro(CDS(cuff),replicates=T)

plot(csDendro(sigGenes,replicates=T), main="Sig Genes")
plot(csDendro(sigTSS,replicates=T), main="Sig TSS")
plot(csDendro(sigCDS,replicates=T), main="Sig CDS")
plot(csDendro(sigIsos,replicates=T), main="Sig Isoforms")

###################################################
### code chunk number 16: global_plots_4
###################################################

ma_g<-MAplot(genes(cuff),over,under) + ggtitle("genes")
ma_g
ma_i<-MAplot(isoforms(cuff),over,under,useCount=T) + ggtitle("isoforms")
ma_i
ma_c<-MAplot(CDS(cuff),over,under,useCount=T) + ggtitle("CDS")
ma_c
ma_t<-MAplot(TSS(cuff),over, under,useCount=T) + ggtitle("TSS")
ma_t

###################################################
### Code Chunk 17:  volcanos
###################################################

vv<-csVolcano(genes(cuff),over, under) + ggtitle("genes")
vv
vy<-csVolcano(TSS(cuff),over, under) + ggtitle("TSS")
vy
vu<-csVolcano(CDS(cuff),over, under) + ggtitle("CDS")
vu
vu<-csVolcano(isoforms(cuff),over, under) + ggtitle("Isoforms")
vu

vol_g<-csVolcanoMatrix(genes(cuff)) + ggtitle("genes")
 vol_g
 vol_i<-csVolcanoMatrix(isoforms(cuff)) + ggtitle("isoforms")
 vol_i
 vol_c<-csVolcanoMatrix(CDS(cuff)) + ggtitle("CDS")
 vol_c
 vol_t<-csVolcanoMatrix(TSS(cuff)) + ggtitle("TSS")
 vol_t

vol_sigT<-csVolcano(sigIsos,over,under, features=T) + ggtitle("Sig Isoforms")
vol_sigT
vol_sigC<-csVolcano(sigCDS,over,under, features=T) + ggtitle("Sig CDS")
vol_sigC
vol_sigI<-csVolcano(sigTSS,over,under, features=T) + ggtitle("Sig TSS_group")
vol_sigI
vol_sigG<-csVolcano(sigGenes,over,under, features=T) + ggtitle("Sig Genes")
vol_sigG
vol_sigmyG<-csVolcano(myGenes,over,under, features=T) + ggtitle("Genes of Interest")
vol_sigmyG

##############################################################
### code chunk number 18: create_geneset/genes_of_interest
###############################################################

gene.counts<-count(genes(cuff))
head(gene.counts)

gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)

###################################################

# the 'geneListString' variable is set in the settings file
# The environment variables are initialized at the begining of this script

myGeneIds = unlist(strsplit(geneListString, " "))
myGeneIds

myGeneSet<-getGenes(cuff,myGeneIds)
#myGeneSet<-getGenes(cuff,myGeneIds), sampleIdList=c(over,under))
myGeneSet

myGenes<-getGenes(cuff,geneIdList=myGeneIds)
myGenes

gene.counts<-count(genes(myGenes))
head(gene.counts)
gene.featurenames<-featureNames(genes(myGenes))
head(gene.featurenames)

GeneSet.counts<-count(genes(myGeneSet))
head(GeneSet.counts)
GeneSet.featurenames<-featureNames(genes(myGeneSet))
head(GeneSet.featurenames)

###################################################
### code chunk number 19: geneset_plots
###################################################

exP.iso<-expressionPlot(myGenes, replicates=F,logMode=T,showErrorbars =T) + ggtitle("Genes of interest")
exP.iso
exp.cds2<-expressionPlot(CDS(myGenes), replicates=T, logMode=T) + ggtitle("CDS of Genes of Interest")
exp.cds2
exp.tss2<-expressionPlot(TSS(myGenes), replicates=T,logMode=T,showErrorbars =F) + ggtitle("TSS of Genes of Interest")
exp.tss2
exp.iso2<-expressionPlot(isoforms(myGenes), replicates=T, logMode=T) + ggtitle("Isoforms of Genes of Interest")
exp.iso2

exP.iso<-expressionPlot(myGeneSet, replicates=F,logMode=T,showErrorbars =T) + ggtitle("Genes of interest")
exP.iso
exp.cds2<-expressionPlot(CDS(myGeneSet), replicates=T, logMode=T) + ggtitle("CDS of Genes of Interest")
exp.cds2
exp.tss2<-expressionPlot(TSS(myGeneSet), replicates=T,logMode=T,showErrorbars =F) + ggtitle("TSS of Genes of Interest")
exp.tss2
exp.iso2<-expressionPlot(isoforms(myGeneSet), replicates=T, logMode=T) + ggtitle("Isoforms of Genes of Interest")
exp.iso2

#####################################################
### code chunk number 20: genes_of_interest_heatmaps
#####################################################

hmiso.rep<-csHeatmap(sigIsos,cluster='both',replicates=T,labRow=F,labCol=T, logMode = T ) + ggtitle("Sig Isoforms")
print(hmiso.rep)
hmiso<-csHeatmap(sigIsos,cluster='both',replicates=F,labRow=F,labCol=T, logMode = T ) + ggtitle("Sig Isoforms")
print(hmiso)

hmtss.rep<-csHeatmap(sigTSS,cluster='both',replicates=T,labRow = F,labCol=T, logMode = T) + ggtitle("Sig TSS_Group")
print(hmtss.rep)
hmtss<-csHeatmap(sigTSS, cluster="row",labRow=F, replicates=F) + ggtitle("sigTSS")
hmtss

hmcds.rep<-csHeatmap(sigCDS,cluster='both',replicates=T,labRow = F,labCol=T, logMode = T) + ggtitle("Sig CDS")
print(hmcds.rep)
hmcds<-csHeatmap(sigCDS,cluster='both',replicates=F,labRow = F,labCol=T, logMode = T) + ggtitle("Sig CDS")
print(hmcds)

hmsigG<-csHeatmap(sigGenes,cluster='both', labRow=F,replicates=F) + ggtitle("Sig Genes")
hmsigG
h.repsigG<-csHeatmap(sigGenes,cluster='both',replicates=T, labRow=F) + ggtitle("Sig Genes")
h.repsigG
h.cs1<-csHeatmap(isoforms(sigGenes), cluster=over, labRow=F, replicates=T) + ggtitle("isoforms(sigGenes)")
h.cs1
h.cs2<-csHeatmap(TSS(sigGenes), cluster=over, replicates=T) + ggtitle("TSS(sigGenes)")
h.cs2
h.cs3<-csHeatmap(CDS(sigGenes), cluster=over,labRow=F, replicates=T) + ggtitle("CDS(sigGenes)")
h.cs3

hg<-csHeatmap(myGenes,cluster="both", replicates=T, labRow=T) + ggtitle("Genes of Interest")
hg
ihiso<-csHeatmap(isoforms(myGenes),cluster='both',labRow=T) + ggtitle("Isoforms of interest")
ihiso
thtss<-csHeatmap(TSS(myGenes),cluster='both',labRow=T) + ggtitle("TSS_Groups of Interest")
thtss
ihcds<-csHeatmap(CDS(myGenes),cluster='both',labRow=T) + ggtitle("CDS_Groups of Interest")
ihcds

hg<-csHeatmap(genes(cuff),cluster="both", replicates=T, labRow=T) + ggtitle("All Genes Features")
hg
ihiso<-csHeatmap(isoforms(cuff),cluster='both',labRow=T) + ggtitle("All Isoforms Features")
ihiso
thtss<-csHeatmap(TSS(cuff),cluster='both',labRow=T) + ggtitle("All TSS_Groups Features")
thtss
ihcds<-csHeatmap(CDS(cuff),cluster='both',labRow=T) + ggtitle("All CDS Features")
ihcds

##################################################
### code chunk number 21: geneset_plots_1.5
###################################################

expbar_i<-expressionBarplot(sigIsos, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed Isofroms")
expbar_i
expbar_t<-expressionBarplot(sigTSS, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed TSS_groups")
expbar_t
expbar_c<-expressionBarplot(sigCDS, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed CDS")
expbar_c
expbar_g3<-expressionBarplot(sigGenes, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed Genes")
expbar_g3
expbar_myg<-expressionBarplot(myGenes) + ggtitle("Genes of Interest")
expbar_myg

####################################
### Code Chunk 22-23:  Gene Loop
####################################

genesArray = c(myGeneIds)
for (i in 1:length(genesArray)) {

    myGeneId = genesArray[i]
    myGeneId
    myGene<-getGene(cuff,geneIdList=myGeneId)
    myGene
    names(myGene)
    myGene[is.na(myGene)]<-c(0)
##########################################################
### Code Chunk 24:  ExpressionPlots/charts - Bar and Pie
###########################################################
 exp.myg<-expressionPlot(myGene, replicates=FALSE)
 print(exp.myg)
 exP.tss.rep<-expressionPlot(isoforms(myGene), replicates=T,logMode=T, drawSummary = F) + ggtitle("Isoform of Genes of Interest")
 exP.tss.rep
 exP.tss<-expressionPlot(TSS(myGene), replicates=T,logMode=T, showErrorbars=T, drawSummary=T) + ggtitle("TSS of Genes of Interest")
 exP.tss
 exP.cds<-expressionPlot(CDS(myGene), replicates=T,logMode=T,showErrorbars =T, drawSummary = T) + ggtitle("CDS of Genes of Interest")
 exP.cds
##########################################################
### Code Chunk 24.5 Pie
###########################################################

gp1<-csPie(myGene,fpkm,level="CDS",replicates=TRUE)
 print(gp2)
 gp2<-csPie(myGene,level="isoforms",replicates=TRUE)
 print(gp2)
 gp3<-csPie(myGene,fpkm,level="TSS",replicates=TRUE)
 print(gp3)
 gp4<-csPie(myGene,level="genes",replicates=TRUE)
 print(gp4)

###################################################
### code chunk number 25: features_2
###################################################
 genetrack<-makeGeneRegionTrack(myGene)
 plotTracks(genetrack)

 trackList<-list()
 myStart<-min(features(myGene)$start)
 myEnd<-max(features(myGene)$end)
 myChr<-unique(features(myGene)$seqnames)

 ideoTrack <- IdeogramTrack(genome = gen_v, chromosome = myChr)
 trackList<-c(trackList,ideoTrack)

 axtrack<-GenomeAxisTrack()
 trackList<-c(trackList,axtrack)

 genetrack<-makeGeneRegionTrack(myGene)
 genetrack

 trackList<-c(trackList,genetrack)

 biomTrack<-BiomartGeneRegionTrack(genome=gen_v,chromosome=as.character(myChr),
                                   start=myStart,end=myEnd,name="ENSEMBL",showId=T)

 trackList<-c(trackList,biomTrack)


 conservation <- UcscTrack(genome = gen_v, chromosome = myChr,
                           track = "Conservation", table = "phyloP100wayAll",
                           from = myStart-2000, to = myEnd+2000, trackType = "DataTrack",
                           start = "start", end = "end", data = "score",
                           type = "hist", window = "auto", col.histogram = "darkblue",
                           fill.histogram = "darkblue", ylim = c(-3.7, 4),
                           name = "Conservation")

 trackList<-c(trackList,conservation)
 plotTracks(trackList,from=myStart-2000,to=myEnd+2000)
}
###################################################
### code chunk number 26: dist_heat_maps
###################################################

DistHeatG<-csDistHeat(genes(cuff),replicates=F) + ggtitle("Genes (w/o replicates)")
DistHeatG
RepDistHeatG<-csDistHeat(genes(cuff),replicates=T) + ggtitle("Genes (w/ replicates)")
RepDistHeatG

DistHeatI<-csDistHeat(isoforms(cuff),replicates=F) + ggtitle("genes (w/o replicates)")
DistHeatI
RepDistHeatI<-csDistHeat(isoforms(cuff),replicates=T) + ggtitle("Genes (w/ replicates)")
RepDistHeatI

DistHeatT<-csDistHeat(TSS(cuff),replicates=FALSE) + ggtitle("CDS, (w/o replicates)")
DistHeatT
RepDistHeatT<-csDistHeat(TSS(cuff),replicates=T) + ggtitle("TSS, (w/ replicates)")
RepDistHeatT

DistHeatC<-csDistHeat(CDS(cuff),replicates=FALSE) + ggtitle("CDS (w/o replicates)")
DistHeatC
RepDistHeatC<-csDistHeat(CDS(cuff),replicates=T) + ggtitle("CDS (w/o replicates)")
RepDistHeatC

Heat_rel<-csDistHeat(relCDS(cuff)) + ggtitle("relCDS")
Heat_rel
Heat_promoters<-csDistHeat(promoters(cuff)) + ggtitle("Promoters")
Heat_promoters
Heat_splicing<-csDistHeat(splicing(cuff)) + ggtitle("Splicing")
Heat_splicing

sigDistHeatG<-csDistHeat(sigGenes) + ggtitle("Genes-sigGenes")
sigDistHeatG
sigRepDistHeatG<-csDistHeat(sigGenes, replicates=T) + ggtitle("Genes-sigGenes")
sigRepDistHeatG
sigHeat_iso_rep<-csDistHeat(isoforms(sigGenes), replicates=T) + ggtitle("RelCDS-sigGenes")
sigHeat_iso_rep
sigHeat_tss_rep<-csDistHeat(TSS(sigGenes), replicates=T) + ggtitle("Promoters-sigGenes")
sigHeat_tss_rep
sigHeat_cds_rep<-csDistHeat(CDS(sigGenes), replicates=T) + ggtitle("Splicing-sigGenes")
sigHeat_cds_rep

sigHeat_rel<-csDistHeat(relCDS(sigGenes)) + ggtitle("RelCDS-sigGenes")
sigHeat_rel
sigHeat_promoters<-csDistHeat(promoters(sigGenes)) + ggtitle("Promoters-sigGenes")
sigHeat_promoters
sigHeat_splicing<-csDistHeat(splicing(sigGenes)) + ggtitle("Splicing-sigGenes")
sigHeat_splicing
###################################################
### code chunk number 28: dim_reduction_1 PCA
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2",scale=T, replicates=FALSE) + ggtitle("Genes")
genes.PCA
rep_genes.PCA<-PCAplot(genes(cuff),"PC1","PC2",scale=T, replicates=TRUE) + ggtitle("Genes, replicates=T")
rep_genes.PCA

cds.PCA<-PCAplot(CDS(cuff),"PC1","PC2",scale =T, replicates=F) + ggtitle("CDS")
cds.PCA
rep_cds.PCA<-PCAplot(CDS(cuff),"PC1","PC2",scale =T, replicates=TRUE) + ggtitle("CDS, replicates=T")
rep_cds.PCA

tss.PCA<-PCAplot(TSS(cuff),"PC1","PC2",replicates=F,scale =T) + ggtitle("TSS")
tss.PCA
rep_tss.PCA<-PCAplot(TSS(cuff),"PC1","PC2",replicates=TRUE,scale =T) + ggtitle("TSS, replicates=T")
rep_tss.PCA

isos.PCA<-PCAplot(isoforms(cuff),"PC1","PC2",replicates=F,scale =T) + ggtitle("Isoforms")
isos.PCA
rep_isos.PCA<-PCAplot(isoforms(cuff),"PC1","PC2", replicates=T, scale=T, pseudocount = 1) + ggtitle("Isoforms, replicates=T")
rep_isos.PCA

###################################################
### code chunk number 29: dim_reduction_2 MDS
###################################################

genes.MDS.rep<-MDSplot(genes(cuff), replicates=T, logMode=T) + ggtitle("Genes, replicates=T")
print(genes.MDS.rep)
genes.MDS.iso<-MDSplot(isoforms(cuff),logMode=T, replicates=TRUE) + ggtitle("Isoforms, replicates=T")
genes.MDS.iso
genes.MDS.cds<-MDSplot(CDS(cuff),logMode=T, replicates=TRUE) + ggtitle("CDS, replicates=T")
genes.MDS.cds
genes.MDS.tss<-MDSplot(TSS(cuff),logMode=T, replicates=TRUE) + ggtitle("TSS, replicates=T")
genes.MDS.tss

siggenes.MDS<-MDSplot(sigGenes,replicates=TRUE, logMode=T) + ggtitle("SIGGENES, replicates=T")
siggenes.MDS
sig_genes.MDS.iso<-MDSplot(isoforms(sigGenes),logMode=T, replicates=TRUE) + ggtitle("SIG_ISOFORMS, replicates=T")
sig_genes.MDS.iso
sig_genes.MDS.cds<-MDSplot(CDS(sigGenes),logMode=T, replicates=TRUE) + ggtitle("SIG_CDS, replicates=T")
sig_genes.MDS.cds
sig_genes.MDS.tss<-MDSplot(TSS(sigGenes),logMode=T, replicates=TRUE) + ggtitle("SIG_TSS, replicates=T")
sig_genes.MDS.tss


#################################################
### code chunk number 29: geneset_cluster_1
###################################################

 ic<-csCluster(genes(myGenes),k=4)
 head(ic$cluster)
 icp<-csClusterPlot(ic)
 icp

###################################################
### code chunk number 30: specificity_1
###################################################

 myGenes.spec<-csSpecificity(genes(myGenes))
 head(myGenes.spec)

#####################################################
##### code chunk number 31: similar_1
#####################################################

 mySimilar<-findSimilar(genes(cuff),myGene,n=20)
 mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)

######################################################
##### code chunk number 32: similar_plots_1
#####################################################

mySimilar<-findSimilar(cuff,diffgene,n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
print(mySimilar.expression)

###################################################
### code chunk number 33: similar_2
###################################################
 myProfile<-c(500,0,400)
 mySimilar2<-findSimilar(CDS(cuff),myProfile,n=10)
 mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
 mySimilar2.expression
 print(mySimilar2.expression)
###################################################
### code chunk number 34: similar_plots_2
###################################################
 mySimilar2.expression<-expressionPlot(mySimilar2,logMode=F,showErrorbars=F)
 mySimilar2.expression
 print(mySimilar2.expression)

###################################################
#gene.fpkm<-fpkm(genes(cuff))
#head(gene.fpkm)
#gene.repFpkm<-repFpkm(genes(cuff))
#head(gene.repFpkm)
#gene.matrix<-fpkmMatrix(genes(cuff))
#head(gene.matrix)
#head(fpkmMatrix(isoforms(cuff))); #head(fpkmMatrix(CDS(cuff)))
#head(fpkmMatrix(TSS(cuff))); #head(fpkmMatrix(sigGenes))
#head(fpkmMatrix(sigIsos)); #head(fpkmMatrix(sigCDS))
#head(fpkmMatrix(sigTSS)); #head(repFpkmMatrix(isoforms(cuff)))
#head(repFpkmMatrix(TSS(cuff))); #head(repFpkmMatrix(CDS(cuff)))
#head(repFpkmMatrix(sigGenes)); #head(repFpkmMatrix(CDS(sigGenes)))
#head(repFpkmMatrix(isoforms(sigGenes))); #head(repFpkmMatrix(TSS(sigGenes)))

###################################################
### code chunk number 9: countMatrix
###################################################
#head(countMatrix(genes(cuff))); #head(countMatrix(TSS(sigGenes)))
#head(countMatrix(isoforms(sigGenes))); #head(countMatrix(CDS(sigGenes)))
#head(countMatrix(sigGenes)); #head(countMatrix(sigIsos))
#head(countMatrix(sigTSS)); #head(countMatrix(sigCDS))


end<-dbDisconnect(cuff@DB)
dev.off(file.path(Rplots))