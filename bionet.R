###########################################################################
# Subnetwork built using bum model of optimal and largest component of connected modules from the human PPI network (interactome)
# Zazil E., Cory C., Todd R.
# Riley lab -UMASS Boston Biology
############################################################################
source("http://bioconductor.org/biocLite.R")
library(BioNet)
library(DLBCL)
## library(RCytoscape)
data(interactome)

options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number, prefixed by #
## options(echo=TRUE)

# We need to keep track of versions used, especially packages
sessionInfo()

####################################################################################
# Evaluate command line arguments ##################################################
####################################################################################

# dir = "~/Documents/UMB/Riley_Lab/HiSeq_data/jose/"
## args = commandArgs(trailingOnly = TRUE)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

geneAlphaValue = as.numeric(Sys.getenv("BIONET_GENE_ALPHA_VALUE"))
subnetAlphaValue = as.numeric(Sys.getenv("BIONET_SUBNET_ALPHA_VALUE"))
numberNetworks = as.numeric(Sys.getenv("BIONET_NUMBER_NETWORKS"))
print(paste("geneAlphaValue = ", geneAlphaValue, sep=""))
print(paste("subnetAlphaValue = ", subnetAlphaValue, sep=""))
print(paste("numberNetworks = ", numberNetworks, sep=""))

####################################################################################
# Functions ########################################################################
####################################################################################

head2 = function(object)
{
    if (is.null(object)) {
        return("NULL")
    }
    return(head(object))
}

reformatSym = function(cuffdiff.table)
{
  # reformat Symbol column to match node labels in the graph
  cuffdiff.table$gene = as.character(cuffdiff.table$gene)
  cuffdiff.table$gene = as.character(cuffdiff.table$gene)
  for (node1 in interactome@nodes) {
    sym = interactome@nodeData@data[node1][[1]]$geneSymbol
    if (!is.null(cuffdiff.table["gene"][cuffdiff.table["gene"] == sym])) {
      cuffdiff.table["gene"][cuffdiff.table["gene"] == sym] = node1
      ## writeLines(c("\nsym=", sym, "; node1=",node1))
    }
  }
  return(cuffdiff.table)
}

createModule = function(aCuffdiffTable)
{
  writeLines(c("\ndim(aCuffdiffTable): ", dim(aCuffdiffTable)))
  writeLines(c("\nlength(aCuffdiffTable$q_value): ", length(aCuffdiffTable$q_value)))

  ## aCuffdiffTable$q_value = p.adjust(aCuffdiffTable$q_value, method = "bonferroni", n = nrow(aCuffdiffTable))
  qval = aCuffdiffTable$q_value # get p values from differential expression data
  ## qval = qval + (1 * 10 ^ -300) #applying correction factor
  names(qval) = aCuffdiffTable$gene #name qvalues with gene symbols

  ## exit

  writeLines(c("\ninteractome:", as.character(summary(interactome))))
  ## interactome

  subnet = subNetwork(aCuffdiffTable$gene, interactome) #Get subnetwork from interactome
  ## writeLines(c("\nsubnet:", head2(subnet)))
  writeLines(c("\nsubnet:", as.character(summary(subnet))))
  ## subnet

  subnet = largestComp(subnet) #include connected nodes only
  ## writeLines(c("\nlargestComp(subnet):", head2(subnet)))
  writeLines(c("\nlargestComp(subnet):", as.character(summary(subnet))))
  ## subnet

  subnet = rmSelfLoops(subnet) #remove self loops
  ## writeLines(c("\nrmSelfLoops(subnet):", head2(subnet)))
  writeLines(c("\nrmSelfLoops(largestComp(subnet)):", as.character(summary(subnet))))
  ## subnet

  bum = fitBumModel(qval, plot = FALSE) # fit Beta-uniform mixture model
  scores = scoreNodes(subnet, bum, fdr = subnetAlphaValue) # score nodes FDR .001 or which?
  logFC = as.numeric(aCuffdiffTable$log2.fold_change.) # old cuffdiff *.diff format
  names(logFC) = aCuffdiffTable$gene
  module = runFastHeinz(subnet, scores) #create optimal module based on scores
  FDR = scanFDR(bum, fdr=c(subnetAlphaValue, subnetAlphaValue/10), labels=names(qval))

  myList = list(module, qval, scores, logFC, FDR)
  return(myList)
}

saveModule = function(module, pval, scores, logFC, dir_out, filename)
{
  #options(bitmapType ='cairo')
  dev.set()
  pathfile = paste(dir_out, filename, sep= "/")
  outfile = paste(pathfile, ".png", sep = "")
  png(file=outfile)
  plotModule(module, scores = scores, diff.expr = logFC)
  plot3dModule(module, diff.or.scores=logFC, red=c("positive"))
  saveNetwork(module, name=filename, file=outfile, type=c("table"))
  saveNetwork(module, name=filename, file=outfile, type=c("XGMML"))
  saveNetwork(module, name=filename, file=outfile, type=c("sif"))
  dev.off()
}



####################################################################################
####################################################################################
####################################################################################


tablePathName = paste(inDir, "/", inFile, sep = "")
writeLines(c("\ntablePathName:", tablePathName))

## cuffdiff.table = read.delim(paste(inDir, inFile, sep = "/"))
cuffdiff.table = read.table(tablePathName, header=TRUE)

## writeLines(c("\ncuffdiff.table:", head2(cuffdiff.table)))
writeLines("\ncuffdiff.table:")
cuffdiff.table[0:10, ]

# reformat the 'Symbol' column so that genes in the data can be matched with nodes in the network
cuffdiff.table = reformatSym(cuffdiff.table)

writeLines(c("\ndim(cuffdiff.table): ", dim(cuffdiff.table))) # cuffdiff.table with all pvalues for differential expression
writeLines(c("\nmax(cuffdiff.table$log2.fold_change.): ", max(cuffdiff.table$log2.fold_change.)))
writeLines(c("\nmin(cuffdiff.table$log2.fold_change.): ", min(cuffdiff.table$log2.fold_change.)))

min(cuffdiff.table$q_value)
max(cuffdiff.table$q_value)

SignificantTable=cuffdiff.table[cuffdiff.table$q_value < geneAlphaValue,]
Significant = createModule(SignificantTable)
Significant[[5]]
saveModule(module=Significant[[1]], pval=Significant[[2]], scores=Significant[[3]], logFC=Significant[[4]], dir_out=outDir, filename="Significant")
dim(Significant)

## UpRegTable = cuffdiff.table[cuffdiff.table$log2.fold_change. > 0,]
UpRegTable = SignificantTable[SignificantTable$log2.fold_change. > 0, ]
#UpRegTable = UpRegTable[UpRegTable$q_value < 1e-200,]
dim(UpRegTable)

writeLines(c("\ndim(UpRegTable): ", dim(UpRegTable))) # UpRegTable with all pvalues for differential expression
writeLines(c("\nmax(UpRegTable$log2.fold_change.): ", max(UpRegTable$log2.fold_change.)))
writeLines(c("\nmin(UpRegTable$log2.fold_change.): ", min(UpRegTable$log2.fold_change.)))

UpRegModule = createModule(UpRegTable)
saveModule(module=UpRegModule[[1]], pval=UpRegModule[[2]], scores=UpRegModule[[3]], logFC=UpRegModule[[4]], dir_out=outDir, filename="UpRegModule")

DownRegTable = SignificantTable[SignificantTable$log2.fold_change. < 0, ]
dim(DownRegTable)
## DownRegTable = cuffdiff.table[cuffdiff.table$log2.fold_change. < 0, ]
## DownRegTable = DownRegTable[DownRegTable$q_value < 1e-200,]
## DownRegTable = DownRegTable[DownRegTable$q_value < 1e-2,]

writeLines(c("\ndim(DownRegTable): ", dim(DownRegTable))) # DownRegTable with all pvalues for differential expression
writeLines(c("\nmax(DownRegTable$log2.fold_change.): ", max(DownRegTable$log2.fold_change.)))
writeLines(c("\nmin(DownRegTable$log2.fold_change.): ", min(DownRegTable$log2.fold_change.)))

DownRegModule = createModule(DownRegTable)
saveModule(module=DownRegModule[[1]], pval=DownRegModule[[2]], scores=DownRegModule[[3]], logFC=DownRegModule[[4]], dir_out=outDir, filename="DownRegModule")




################### RCytoscape #######################################################
# RCytoscape requires a port to be listening using rpc protocol. Default is port 9000.
# Need to see what happens running the R script in the cluster and Cytoscape in my laptop.
# Most likely host will need to be changed to my IP address.
#######################################################################################
# 3,104 red (positive) nodes    # 6 green (negative) nodes    #0 nodes are white
#sum(scores[nodes(module)] > 0)   #sum(scores[nodes(module)] < 0) #sum(scores[nodes(module)]==0)

#g<-new('graphNEL', edgemode='directed')
#cw<-new.CytoscapeWindow(title="Jose's data", graph=g, overwriteWindow=TRUE, host = "localhost", rpcPort = 9000, )

#g<-initNodeAttribute (module, attribute.name='geneID', attribute.type='numeric', default.value=0)
#g<-initNodeAttribute (g, 'geneSymbol', 'char', 'x')
#g<-initNodeAttribute (g, 'score', 'numeric', scores)
#g<-initNodeAttribute(g,'lfc','numeric',logFC)

##g<-initEdgeAttribute(g,"weight",'integer',0)
##g<-initEdgeAttribute (g, "similarity", 'numeric', 0)
#nodeData(g,cuffdiff.table$gene,"geneID")<-cuffdiff.table$gene_id
#nodeData(g,cuffdiff.table$gene,"geneSymbol")<-cuffdiff.table$gene
#nodeData(g,cuffdiff.table$gene,"score")<-scores
#nodeData(g,cuffdiff.table$gene,"lfc")<-logFC
#cw=setGraph(cw,g)

#set node colors based on logFC
#green when logFC is negative, red is positive, white is zero and shades are in between
#setNodeColorRule(cw,"lfc", c (-3.0, 0.0, 3.0), c('#00AA00', '#00FF00', '#FFFFFF', '#FF0000', '#AA0000'), mode='interpolate')
#displayGraph(cw)
#redraw(cw)
#layoutNetwork(cw, hlp[18]) #default fruchterman.reingold layout is 18

# Apply values to some of the properties and plot the layout
#cy <- CytoscapeConnection()
#hlp <-getLayoutNames(cy)

#writeLines(c(noa.names(getGraph(cw)))) #what data attributes are defined? noa=node attributes
##writeLines(c(noa(getGraph(cw),'geneID')))
##writeLines(c(noa(getGraph(cw),'geneSymbol')))
##writeLines(c(noa(getGraph(cw),'score')))
##writeLines(c(noa(getGraph(cw),'lfc')))

#saveNetwork(g, "Jose_RCytoscape", 'XGMML')  #doesn't work, so I just saved it directly from Cytoscape
############################################################

