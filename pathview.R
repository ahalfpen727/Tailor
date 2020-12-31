source("http://bioconductor.org/biocLite.R")

options(bitmapType='cairo')
options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number, prefixed by #
## options(echo=TRUE)

sessionInfo() # Keep track of versions used, especially packages


alphaDown = as.numeric(Sys.getenv("PATHVIEW_ALPHA_DOWN"))
alphaUp = as.numeric(Sys.getenv("PATHVIEW_ALPHA_UP"))
species = Sys.getenv("PATHVIEW_SPECIES")
numberPathwaysDown = as.numeric(Sys.getenv("PATHVIEW_NUMBER_PATHWAYS_DOWN"))
numberPathwaysUp = as.numeric(Sys.getenv("PATHVIEW_NUMBER_PATHWAYS_UP"))
pathwaysList = Sys.getenv("PATHVIEW_PATHWAYS_LIST")

## test = Sys.getenv("PATH")
## print(paste("test = ", test, sep=""))
## test = Sys.getenv("QUEUE")
## print(paste("test = ", test, sep=""))

print(paste("alphaDown = ", alphaDown, sep=""))
print(paste("alphaUp = ", alphaUp, sep=""))
print(paste("species = ", species, sep=""))
print(paste("numberPathwaysDown = ", numberPathwaysDown, sep=""))
print(paste("numberPathwaysUp = ", numberPathwaysUp, sep=""))
print(paste("pathwaysList = ", pathwaysList, sep=""))


####################################################################################
# Evaluate command line arguments ##################################################
####################################################################################
# dir = "~/Documents/UMB/Riley_Lab/HiSeq_data/jose/"
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

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

####################################################################################
####################################################################################
####################################################################################


## cuffdiff.table = read.xlsx2(paste(inDir, inFile, sep = ""), startRow = 2)
tablePathName = paste(inDir, "/", inFile, sep = "")
writeLines(c("\ntablePathName:", tablePathName))

cuffdiff.table = read.table(tablePathName, header=TRUE)

## writeLines("\ncuffdiff.table:", head2(cuffdiff.table))
writeLines("\ncuffdiff.table:")
cuffdiff.table[0:10, ]

## biocLite("pathview")
## biocLite("gage")

outDir
setwd(outDir)
getwd()

require(pathview)
require(gage)

simpleCap <- function(x) {
  x <- tolower(x)
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

## cuffdiff.table = read.delim(file = "gene_exp.diff", sep = "\t")

## Cufflinks is a very popular RNA-Seq analysis tool.  It is developed by the same group as TopHat.  It is
## implemented independent of Bioconductor.  However, we can read its differential expression analysis
## results into R easily.  The result file is named *_exp_diff
## Notice that the gene symbols need to be converted to Entrez Gene IDs, which are used in
## KEGG pathways (for many research species) and GO gene sets.

# notice the column name special character changes. The column used to be
# cuffdiff.table$log2.fold_change. for older versions of Cufflinks.
# cuff.fc = cuffdiff.table$log2.FPKMy.FPKMx
cuff.fc = cuffdiff.table$log2.fold_change.

writeLines(c("\ncuff.fc: ", head2(cuff.fc)))

gnames = cuffdiff.table$gene
sel = gnames != "-"
head(sel)
## writeLines("sel:")
## sel[0:10]

gnames = as.character(gnames[sel])
head(gnames)

## writeLines("gnames:")
## gnames[0:10]
## head(gnames)

cuff.fc = cuff.fc[sel]
head(cuff.fc)
## cuff.fc[0:10]

names(cuff.fc) = gnames
## gnames[0:10]

# Notice that the gene sybmols need to be converted to Entrez Gene IDs, which are used in
                                        # KEGG pathways (many research species) and GO gene sets.
data(korg)
## labels = colnames(korg)
## labels[3] = "commonName"
## labels[1] = kegg.code
## colnames(korg) = labels
head(korg)
dim(korg)
colnames(korg)

kegg.code = korg[which(korg[,"scientific.name"] == species),][["kegg.code"]]
kegg.code

# ADDED LINES

# The following line may serve as a substitute for species not in the conversational lexicon and lacking common names (i.e. bacteria)
# kegg.code = korg[which(korg[,"kegg.code"] == species),][[2]]; kegg.code
# kegg.code = korg[which(korg[,"scientific.name"] == species),][[1]]; kegg.code
# kegg.code = korg[which(korg[,"common.name"] == species),][[1]]

org.code = substr(simpleCap(kegg.code),1, 2)
org.code

gnames.entrez = pathview::id2eg(gnames, category = "symbol", org = org.code)
head(gnames.entrez)

#Modified Line
#gnames.entrez = pathview::id2eg(gnames, category = "symbol", org = kegg.code)

sel2 = gnames.entrez[,2]>""
head(sel2)

cuff.fc = cuff.fc[sel2]
names(cuff.fc) = gnames.entrez[sel2,2]
writeLines(c("\nrange(cuff.fc): ", range(cuff.fc)))

# remove the -Inf and Inf values, which block the downstream analysis
cuff.fc[cuff.fc > 10] = 10
cuff.fc[cuff.fc < -10] = -10

exp.fc = cuff.fc[ !is.na(cuff.fc) ]

myPathviewIds = unlist(strsplit(pathwaysList, " "))

if (length(myPathviewIds) > 0) {
    pv.out.list <- sapply(myPathviewIds,
                          function(pid) pathview(gene.data =  exp.fc, pathway.id = pid,  species = kegg.code, out.suffix = "cuffdiff", kegg.dir = ".", kegg.native=TRUE))

    pv.out.list
}
else {

    ## Next, we use GAGE for pathway analysis, and Pathview for visualization.  Notice that
    ## this step (the same code) is identical for DESeq, edgeR, Limma and Cufflinks workflows.

                                        # In pathview function, please be conscious of the argument gene.idtype. Default is gene.idtype="entrez",i.e. Entrez
                                        # Gene, which are the primary KEGG gene ID for many common model organisms. For other species, gene.idtype
                                        # should be set to "KEGG" since KEGG uses other types of gene IDs.

    ## data(kegg.gs)
    kg.species<- kegg.gsets(species)
    kegg.gs<- kg.species$kg.sets[kg.species$sigmet.idx]

                                        #let's get the annotation files for mouse and convert the gene set to gene symbol format
    writeLines(c("\nlapply(kegg.gs[1:3], head):", lapply(kegg.gs[1:3], head)))
    writeLines(c("\nlength(kegg.gs): ", length(kegg.gs)))
    writeLines(c("\nexp.fc:", head2(exp.fc)))
    writeLines(c("\nstr(exp.fc):", str(exp.fc)))

    fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)

    writeLines(c("\nfc.kegg.p$greater:", head2(fc.kegg.p$greater)))

    writeLines(c("\nfc.kegg.p$greater[, \"p.val\"]:", (fc.kegg.p$greater[, "p.val"])))

    upExpressed <- fc.kegg.p$greater[, "p.val"] < alphaUp & !is.na(fc.kegg.p$greater[, "p.val"])
    downExpressed <- fc.kegg.p$less[, "p.val"] < alphaDown & !is.na(fc.kegg.p$less[,"p.val"])

    print(paste("max(upRegulated p.val) = ", max(fc.kegg.p$greater[, "p.val"]),  sep=""))
    print(paste("min(upRegulated p.val) = ", min(fc.kegg.p$greater[, "p.val"]),  sep=""))

    print(paste("max(downRegulated p.val) = ", max(fc.kegg.p$less[, "p.val"]), sep=""))
    print(paste("min(downRegulated p.val) = ", min(fc.kegg.p$less[, "p.val"]),  sep=""))

    path.ids.upExpr <- rownames(fc.kegg.p$greater)[upExpressed]
    path.ids2.upExpr <- substr(path.ids.upExpr, 1, 8)
    writeLines(c("\npath.ids.upExpr:", head2(path.ids.upExpr)))
    writeLines(c("\nlength(path.ids.upExpr):", length(path.ids.upExpr)))

    path.ids.downExpr <- rownames(fc.kegg.p$less)[downExpressed]
    path.ids2.downExpr <- substr(path.ids.downExpr, 1, 8)
    writeLines(c("\npath.ids.downExpr:", head2(path.ids.downExpr)))
    writeLines(c("\nlength(path.ids.downExpr):", length(path.ids.downExpr)))

    ## path.ids.upOrDownExpr <- substr(c(path.ids, path.ids.l), 1, 8)
    ## writeLines(c("\npath.ids.upOrDownExpr:", head2(path.ids.upOrDownExpr)))
    ## writeLines(c("\nlength(path.ids.upOrDownExpr):", length(path.ids.upOrDownExpr)))


                                        # view first N pathways as demo
                                        # species = "hsa" is for human
    topN = min(numberPathwaysDown, length(path.ids.upExpr))
    topN
    if (topN > 0) {
        pv.out.list <- sapply(path.ids2.upExpr[1:topN],
                              function(pid) pathview(gene.data =  exp.fc, pathway.id = pid,  species = kegg.code, out.suffix = "cuffdiff", kegg.dir = ".", kegg.native=TRUE))

        pv.out.list
    }

    topN = min(numberPathwaysUp, length(path.ids.downExpr))
    topN
    if (topN > 0) {
        pv.out.list <- sapply(path.ids2.downExpr[1:topN],
                              function(pid) pathview(gene.data =  exp.fc, pathway.id = pid,  species = kegg.code, out.suffix = "cuffdiff", kegg.dir = ".", kegg.native=TRUE))

        pv.out.list
    }
    ## topN = min(numberPathways, length(path.ids.upOrDownExpr))
    ## pv.out.list <- sapply(path.ids.upOrDownExpr[1:topN],
    ##                       function(pid) pathview(gene.data =  exp.fc, pathway.id = pid,  species = kegg.code, out.suffix = "cuffdiff", kegg.dir = outDir, kegg.native=TRUE))

## pv.out.list
}

## Here we used GAGE to infer the significant pathways.  But we are not limited to these pathways.  We can
## use Pathview to visualize RNA-Seq data (exp.fc here) on all interesting pathways directly.  We may also do GO
## and other types of gene set analysis as described in native workflow above.

