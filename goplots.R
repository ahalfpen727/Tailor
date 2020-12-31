###############################################################################
# GOplots for Significantly differentially expressed features: method 1
############################################################################
# Libraries
###########################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO");biocLite("org.Hs.eg.db");biocLite("ReactomePA")
#biocLite("pathview");biocLite("reactome.db");biocLite("clusterProfiler")
#biocLite("GSEABase");biocLite("DOSE")
#biocLite("GO.db");biocLite("KEGG.db")

library(cummeRbund);library(biomaRt);library(reactome.db);library(ReactomePA)
library(org.Mm.eg.db);library(org.Hs.eg.db)
library(topGO);library(GSEABase);library(clusterProfiler)
library(GO.db);library(KEGG.db);library(DOSE); library(pathview)
###############################################################################
# GOplots for Significantly differentially expressed features: method 1
############################################################################
## source('openGraphSaveGraph.r');
options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number
# prefixed by #

### original  R code from Drew'
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

###################################################
### code chunk number 1: init, 2: loadLib, 3: Read
###################################################
options(width=65)
library(cummeRbund)
refgtf = Sys.getenv("CUMMERBUND_INPUT_GTF")
genome_path = Sys.getenv("GENOME_PATH")
gen_v = Sys.getenv("R_GENOME")
cuff<-readCufflinks(gtfFile=refgtf,genome=genome_path,rebuild=T)
cuff
###################################################
### for testing
###################################################
# default
#cuff<-readCufflinks(gtfFile="GRCh38_genes.gtf",genome="genome.fa",rebuild=F)
#cuff
#over="LUTS"
#under="CTRL"

###############################################################################
# Features, Counts, and all inclusive tables
# get gene symbols for CDS TSS Isoforms features
################################################################################

runInfo(cuff);
repdata<-replicates(cuff)
reps<-repdata$rep_name
mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
sigGenes<-getGenes(cuff, mySigGenes)
length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)
sig_genes_annot<-sigGenes@annotation
head(sig_genes_exp.diff)
head(sig_genes_annot)
sig_genes_exp.diff["gene_id"]<-sig_genes_annot["gene_short_name"]
head(sig_genes_exp.diff)

genes_exp.diff<-diffData(isoforms(cuff))
dim(genes_exp.diff)
g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(g.rep.matrix)
gene.xloc.matrix<-featureNames(isoforms(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list@isoforms)
dim(gene_annotation_data)
gene.rep.matrix<-cbind2(gene_annotation_data, g.rep.matrix)
gene_exp.diff<-cbind2(gene_annotation_data, genes_exp.diff)
genes_exp.diff<-as.data.frame(cbind(gene_annotation_data, logFC=gene_exp.diff$log2_fold_change,pvalue=gene_exp.diff$p_value, qvalue=gene_exp.diff$q_value))

gene_exp.diff<-genes_exp.diff[unique(genes_exp.diff$gene_short_name),]
repdata<-replicates(cuff)
reps<-repdata$rep_name
samples<-repdata$sample_name
colnames(g.rep.matrix)<-samples
groups<-factor(g.rep.matrix,levels=0:1,labels=c(over, under))
g.over.matrix<-g.rep.matrix[samples == over]
g.under.matrix<-g.rep.matrix[samples == under]

colnames(g.rep.matrix)<-reps
rownames(g.rep.matrix)<-gene_annotation_data$tracking_id
genes.reps.df<-cbind(gene_annotation_data$gene_short_name, g.rep.matrix)
genes.reps.df<-gene.rep.matrix[unique(gene.rep.matrix$gene_short_name),]
genes<-genes.reps.df$gene_short_name

foldchange<-sig_genes_exp.diff[,"log2_fold_change"]
under.inf<-which(foldchange == "-Inf")
over.inf = which(foldchange == "Inf")
over.= which(foldchange > 0, foldchange == "Inf", foldchange != "-Inf")
under. = which(foldchange < 0,  foldchange == "-Inf")

NOexp.inUNDER<-as.data.frame(sig_genes_exp.diff[over.inf,])
NOexp.inOVER<-as.data.frame(sig_genes_exp.diff[under.inf,])

LOexp.inUNDER<-as.data.frame(sig_genes_exp.diff[over.,])
LOexp.inOVER<-as.data.frame(sig_genes_exp.diff[under.,])

QvalnoOVER<-NOexp.inOVER[,"q_value"]
names(QvalnoOVER)<-NOexp.inOVER[,"gene_id"]

QvalnoUNDER<-NOexp.inUNDER[,"q_value"]
names(QvalnoUNDER)<-NOexp.inUNDER[,"gene_id"]

OVERupexp<-as.data.frame(sig_genes_exp.diff[over.,])
QvalOVERhi<-OVERupexp[,"q_value"]
names(QvalOVERhi)<-OVERupexp[,"gene_id"]

UNDERupexp<-as.data.frame(sig_genes_exp.diff[under.,])
QvalUNDERhi<-UNDERupexp[,"q_value"]
names(QvalUNDERhi)<-UNDERupexp[,"gene_id"]

#QvalOVERhi #QvalUNDERhi
ENTREZQvalUNDERhi<-mapIds(x = org.Hs.eg.db,
                          keys = names(QvalUNDERhi),
                          column = "ENTREZID",
                          keytype = "SYMBOL",
			   multiVals="first")

ENTREZQvalOVERhi<-mapIds(x = org.Hs.eg.db,
                          keys = names(QvalOVERhi),
                          column = "ENTREZID",
                          keytype = "SYMBOL",
			  multiVals="first")

ENTREZsiggenes<-mapIds(x = org.Hs.eg.db,
                          keys = sig_genes_exp.diff$gene_id,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
			  multiVals="first")

under.hilog<-sig_genes_exp.diff[,"log2_fold_change"]
names(under.hilog)<-ENTREZsiggenes
under.hiqval<-sig_genes_exp.diff[,"q_value"]
names(under.hiqval)<-ENTREZQvalUNDERhi
under.hipval<-(sig_genes_exp.diff[,"p_value"])
names(under.hipval) <- ENTREZQvalUNDERhi

over.hilog<-c(sig_genes_exp.diff[,"log2_fold_change"])
names(over.hilog)<-(ENTREZQvalOVERhi)
over.hiqval<-sig_genes_exp.diff[,"q_value"]
names(over.hiqval)<-ENTREZQvalOVERhi
over.hipval<-sig_genes_exp.diff[,"p_value"]
names(over.hipval)<-ENTREZQvalOVERhi

ENTREZQvalUNDER<-mapIds(x = org.Hs.eg.db,
                          keys =NOexp.inUNDER[,"gene_id"],
                          column = "ENTREZID",
                          keytype = "SYMBOL",
							multiVals="first")

ENTREZQvalOVER<-mapIds(x = org.Hs.eg.db,
                          keys =NOexp.inOVER[,"gene_id"],
                          column = "ENTREZID",
                          keytype = "SYMBOL",
						       multiVals="first" )

ENTREZids<-mapIds(x = org.Hs.eg.db,
                          keys = gene.list@isoforms@annotation$gene_short_name,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
					          multiVals="first")

sig_genes_exp<-cbind(ENTREZsiggenes,log2fold=sig_genes_exp.diff[,"log2_fold_change"])
head(sig_genes_exp)
dim(sig_genes_exp)
names(sig_genes_exp) = rownames(sig_genes_exp)

siggenelog<-sig_genes_exp.diff[,"log2_fold_change"]
names(siggenelog)<-ENTREZsiggenes
siggeneqval<-sig_genes_exp.diff[,"q_value"]
q.val<-sig_genes_exp.diff[,"q_value"]
names(siggeneqval)<-q.val
names(siggeneqval)<-ENTREZsiggenes
siggenepval<-sig_genes_exp.diff[,"p_value"]
names(siggenepval)<-ENTREZsiggenes

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
topDiffGenes <- function(allScore) {
return(allScore < 0.01)}

# Significantly Differentially Expressed by the Under group
# i.e Luts over Ctrl --> Ctrl = Under group, Under group highly
# expresses and No expression in the Over group

noOVERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalnoOVER, nodeSize = 5,
	       	  			     annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					     	       ID = "symbol")
under.upBPGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalUNDERhi, nodeSize = 5,
		    			       annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					       	         ID = "symbol")
noOVERmfGOdata <- new("topGOdata",ontology = "MF",allGenes = QvalnoOVER, nodeSize = 5,
	       	  			     annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					     	       ID = "symbol")
under.upMFGOdata <- new("topGOdata",ontology = "MF",allGenes = QvalUNDERhi, nodeSize = algorithm = "elim", statistic = "ks")

over.upBPtKS <- runTest(over.upBPGOdata, algorithm = "classic", statistic = "ks")
over.upBPFisher <- runTest(over.upBPGOdata, algorithm = "classic", statistic = "fisher")
over.upBPtKS.elim <- runTest(over.upBPGOdata, algorithm = "elim", statistic = "ks")



over.upBPtKS <- runTest(over.upBPGOdata, algorithm = "classic", statistic = "ks")
over.upBPFisher <- runTest(over.upBPGOdata, algorithm = "classic", statistic = "fisher")
over.upBPtKS.elim <- runTest(over.upBPGOdata, algorithm = "elim", statistic = "ks")

over.upMFtKS <- runTest(over.upMFGOdata, algorithm = "classic", statistic = "ks")
over.upMFFisher <- runTest(over.upMFGOdata, algorithm = "classic", statistic = "fisher")
over.upMFtKS.elim <- runTest(over.upMFGOdata, algorithm = "elim", statistic = "ks")

over.upCCtKS <- runTest(over.upCCGOdata, algorithm = "classic", statistic = "ks")
over.upCCFisher <- runTest(over.upCCGOdata, algorithm = "classic", statistic = "fisher")
over.upCCtKS.elim <- runTest(over.upCCGOdata, algorithm = "elim", statistic = "ks")

pdf("GOplots_OVER_^.pdf")

showSigOfNodes(over.upMFGOdata, score(over.upMFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(over.upMFGOdata, score(over.upMFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(over.upMFGOdata, score(over.upMFtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(over.upMFGOdata, over.upMFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(over.upMFGOdata, over.upMFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(over.upMFGOdata, over.upMFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(over.upCCGOdata, score(over.upCCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(over.upCCGOdata, score(over.upCCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(over.upCCGOdata, score(over.upCCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(over.upCCGOdata, over.upCCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(over.upCCGOdata, over.upCCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(over.upCCGOdata, over.upCCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(over.upBPGOdata, score(over.upBPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(over.upBPGOdata, score(over.upBPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(over.upBPGOdata, score(over.upBPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(over.upBPGOdata, over.upBPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(over.upBPGOdata, over.upBPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(over.upBPGOdata, over.upBPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noUNDERccGOdata, score(noUNDERccFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noUNDERccGOdata, score(noUNDERcctKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(noUNDERccGOdata, score(noUNDERcctKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(noUNDERccGOdata, noUNDERccFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noUNDERccGOdata, noUNDERcctKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noUNDERccGOdata, noUNDERcctKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noUNDERbpGOdata, score(noUNDERbpFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noUNDERbpGOdata, score(noUNDERbptKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(noUNDERbpGOdata, score(noUNDERbptKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(noUNDERbpGOdata, noUNDERbpFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noUNDERbpGOdata, noUNDERbptKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noUNDERbpGOdata, noUNDERbptKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noUNDERmfGOdata, score(noUNDERmfFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noUNDERmfGOdata, score(noUNDERmftKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(noUNDERmfGOdata, score(noUNDERmftKS.elim), firstSigNodes = 5, useInfo = "def" )
pri