# This GO analysis uses the GOstats package with a list of significantly differentially expressed genes.
# However, this analysis does not use the p-values or q-values of the differentially expressed genes.


source("http://bioconductor.org/biocLite.R")


options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number, prefixed by #
## options(echo=TRUE)

# We need to keep track of versions used, especially packages
sessionInfo()

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
print(args)
print(R_LIBS_USER)
print(R_LIBS)
####################################################################################


diffExprAlpha = as.numeric(Sys.getenv("GO_DIFF_EXPRESSED_ALPHA_VALUE"))
hyperGeoAlpha = as.numeric(Sys.getenv("GO_HYPER_GEO_ALPHA_VALUE"))
print(paste("diffExprAlpha = ", diffExprAlpha, sep=""))
print(paste("hyperGeoAlpha = ", hyperGeoAlpha, sep=""))

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

# read in tables
## library(xlsx)

# increase jvm memory, otherwise it crashes half the time
# options(java.parameters = "-Xmx3072m")

## table = read.xlsx2(paste(inDir, inFile, sep = ""), startRow = 2)
tablePathName = paste(inDir, "/", inFile, sep = "")
writeLines(c("\ntablePathName:", tablePathName))

cuffdiff.table = read.table(tablePathName, header=TRUE)

## writeLines(c("\ncuffdiff.table before filtering:", head2(cuffdiff.table)))
writeLines("\ncuffdiff.table before filtering:")
cuffdiff.table[0:10, ]

# include only significantly differentially expressed genes
## cutoff = 1 * (10 ^ -1)
cuffdiff.table = cuffdiff.table[cuffdiff.table$q_value < diffExprAlpha, ]

## writeLines(c("\ncuffdiff.table after filtering:", head2(cuffdiff.table)))
writeLines("\ncuffdiff.table after filtering:")
cuffdiff.table[0:10, ]

#biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

#biocLite("GOstats")
library(GOstats)
require(pathview)

# Notice that the gene sybmols need to converted to Entrez Gene IDs, which are used in
# KEGG pathways (many research species) and GO gene sets.
gnames.eg = pathview::id2eg(cuffdiff.table$gene, category ="symbol")
sel2 = gnames.eg[,2] > ""
gnames.eg = gnames.eg[sel2]

## writeLines("cuffdiff.table$gene_id:")
## cuffdiff.table$gene_id[0:10]

## writeLines("sel2:")
## sel2[0:10]


entrez_object = org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes = mappedkeys(entrez_object)

##################################################################################################
# Over-representation
##################################################################################################


# Cellular component
params.over.uncon.cc = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
)
hgOver.uncon.cc = hyperGTest(params.over.uncon.cc)
result.over.cc = summary(hgOver.uncon.cc)
result.over.cc$Adj.Pvalue = p.adjust(result.over.cc$Pvalue, method = 'BH')
result.over.cc = cbind(result.over.cc[1:2], result.over.cc["Adj.Pvalue"], result.over.cc[3:7])
writeLines(c("\nresult.over.cc:", head2(result.over.cc)))

# Cellular component
params.over.con.cc = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=T,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
)
hgOver.cc = hyperGTest(params.over.con.cc)
result.over.cc = summary(hgOver.cc)
result.over.cc$Adj.Pvalue = p.adjust(result.over.cc$Pvalue, method = 'BH')
result.over.cc = cbind(result.over.cc[1:2], result.over.cc["Adj.Pvalue"], result.over.cc[3:7])
writeLines(c("\nresult.over.cc:", head2(result.over.cc)))

# Biological process
params.over.uncon.bp = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('BP'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
)
hgOver.uncon.bp = hyperGTest(params.over.uncon.bp)
result.over.bp = summary(hgOver.uncon.bp)
result.over.bp$Adj.Pvalue = p.adjust(result.over.bp$Pvalue, method = 'BH')
result.over.bp = cbind(result.over.bp[1:2], result.over.bp["Adj.Pvalue"], result.over.bp[3:7])
writeLines(c("\nresult.over.bp:", head2(result.over.bp)))

# Biological process
params.over.con.bp = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('BP'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=T,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
)
hgOver.bp = hyperGTest(params.over.con.bp)
result.over.bp = summary(hgOver.bp)
result.over.bp$Adj.Pvalue = p.adjust(result.over.bp$Pvalue, method = 'BH')
result.over.bp = cbind(result.over.bp[1:2], result.over.bp["Adj.Pvalue"], result.over.bp[3:7])
writeLines(c("\nresult.over.bp:", head2(result.over.bp)))

# Molecular function
params.over.uncon.mf = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
)
hgOver.uncon.mf = hyperGTest(params.over.uncon.mf)
result.over.mf = summary(hgOver.uncon.mf)
result.over.mf$Adj.Pvalue = p.adjust(result.over.mf$Pvalue, method = 'BH')
result.over.mf = cbind(result.over.mf[1:2], result.over.mf["Adj.Pvalue"], result.over.mf[3:7])
writeLines(c("\nresult.over.mf:", head2(result.over.mf)))

# Molecular function
params.over.con.mf = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=T,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
)
hgOver.mf = hyperGTest(params.over.con.mf)
result.over.mf = summary(hgOver.mf)
result.over.mf$Adj.Pvalue = p.adjust(result.over.mf$Pvalue, method = 'BH')
result.over.mf = cbind(result.over.mf[1:2], result.over.mf["Adj.Pvalue"], result.over.mf[3:7])
writeLines(c("\nresult.over.mf:", head2(result.over.mf)))

# HTML doc results
htmlReport(hgOver.uncon.bp, paste(outDir, "/go_biologicalProcess_uncon.over.html", sep = ""))
htmlReport(hgOver.uncon.cc, paste(outDir, "/go_cellularComponent_uncon.over.html", sep = ""))
htmlReport(hgOver.uncon.mf, paste(outDir, "/go_molecularFunction_uncon.over.html", sep = ""))
htmlReport(hgOver.bp, paste(outDir, "/go_biologicalProcess_con.over.html", sep = ""))
htmlReport(hgOver.cc, paste(outDir, "/go_cellularComponent_con.over.html", sep = ""))
htmlReport(hgOver.mf, paste(outDir, "/go_molecularFunction_con.over.html", sep = ""))

# txt table of results
write.table(result.over.bp, paste(outDir, "/go_biologicalProcess_uncon.over.txt", sep = ""))
write.table(result.over.cc, paste(outDir, "/go_cellularComponent_uncon.over.txt", sep = ""))
write.table(result.over.mf, paste(outDir, "/go_molecularFunction_uncon.over.txt", sep = ""))
write.table(result.over.bp, paste(outDir, "/go_biologicalProcess_con.over.txt", sep = ""))
write.table(result.over.cc, paste(outDir, "/go_cellularComponent_con.over.txt", sep = ""))
write.table(result.over.mf, paste(outDir, "/go_molecularFunction_con.over.txt", sep = ""))


##################################################################################################
# Under-representation
##################################################################################################


# Cellular component
params.under.uncon.cc = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=F,
                   testDirection='under',
                   annotation="org.Hs.eg.db"
)
hgUnder.uncon.cc = hyperGTest(params.under.uncon.cc)
result.under.cc = summary(hgUnder.uncon.cc)
result.under.cc$Adj.Pvalue = p.adjust(result.under.cc$Pvalue, method = 'BH')
result.under.cc = cbind(result.under.cc[1:2], result.under.cc["Adj.Pvalue"], result.under.cc[3:7])
writeLines(c("\nresult.under.cc:", head2(result.under.cc)))

# Cellular component
params.under.con.cc = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=T,
                   testDirection='under',
                   annotation="org.Hs.eg.db"
)
hgUnder.cc = hyperGTest(params.under.con.cc)
result.under.cc = summary(hgUnder.cc)
result.under.cc$Adj.Pvalue = p.adjust(result.under.cc$Pvalue, method = 'BH')
result.under.cc = cbind(result.under.cc[1:2], result.under.cc["Adj.Pvalue"], result.under.cc[3:7])
writeLines(c("\nresult.under.cc:", head2(result.under.cc)))

# Biological process
params.under.uncon.bp = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('BP'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=F,
                   testDirection='under',
                   annotation="org.Hs.eg.db"
)
hgUnder.uncon.bp = hyperGTest(params.under.uncon.bp)
result.under.bp = summary(hgUnder.uncon.bp)
result.under.bp$Adj.Pvalue = p.adjust(result.under.bp$Pvalue, method = 'BH')
result.under.bp = cbind(result.under.bp[1:2], result.under.bp["Adj.Pvalue"], result.under.bp[3:7])
writeLines(c("\nresult.under.bp:", head2(result.under.bp)))

# Biological process
params.under.con.bp = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('BP'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=T,
                   testDirection='under',
                   annotation="org.Hs.eg.db"
)
hgUnder.bp = hyperGTest(params.under.con.bp)
result.under.bp = summary(hgUnder.bp)
result.under.bp$Adj.Pvalue = p.adjust(result.under.bp$Pvalue, method = 'BH')
result.under.bp = cbind(result.under.bp[1:2], result.under.bp["Adj.Pvalue"], result.under.bp[3:7])
writeLines(c("\nresult.under.bp:", head2(result.under.bp)))

# Molecular function
params.under.uncon.mf = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=F,
                   testDirection='under',
                   annotation="org.Hs.eg.db"
)
hgUnder.uncon.mf = hyperGTest(params.under.uncon.mf)
result.under.mf = summary(hgUnder.uncon.mf)
result.under.mf$Adj.Pvalue = p.adjust(result.under.mf$Pvalue, method = 'BH')
result.under.mf = cbind(result.under.mf[1:2], result.under.mf["Adj.Pvalue"], result.under.mf[3:7])
writeLines(c("\nresult.under.mf:", head2(result.under.mf)))

# Molecular function
params.under.con.mf = new('GOHyperGParams',
                   geneIds=gnames.eg,
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=hyperGeoAlpha,
                   conditional=T,
                   testDirection='under',
                   annotation="org.Hs.eg.db"
)
hgUnder.mf = hyperGTest(params.under.con.mf)
result.under.mf = summary(hgUnder.mf)
result.under.mf$Adj.Pvalue = p.adjust(result.under.mf$Pvalue, method = 'BH')
result.under.mf = cbind(result.under.mf[1:2], result.under.mf["Adj.Pvalue"], result.under.mf[3:7])
writeLines(c("\nresult.under.mf:", head2(result.under.mf)))

# HTML doc results
htmlReport(hgUnder.uncon.bp, paste(outDir, "/go_biologicalProcess_uncon.under.html", sep = ""))
htmlReport(hgUnder.uncon.cc, paste(outDir, "/go_cellularComponent_uncon.under.html", sep = ""))
htmlReport(hgUnder.uncon.mf, paste(outDir, "/go_molecularFunction_uncon.under.html", sep = ""))
htmlReport(hgUnder.bp, paste(outDir, "/go_biologicalProcess_con.under.html", sep = ""))
htmlReport(hgUnder.cc, paste(outDir, "/go_cellularComponent_con.under.html", sep = ""))
htmlReport(hgUnder.mf, paste(outDir, "/go_molecularFunction_con.under.html", sep = ""))

# txt table of results
write.table(result.under.bp, paste(outDir, "/go_biologicalProcess_uncon.under.txt", sep = ""))
write.table(result.under.cc, paste(outDir, "/go_cellularComponent_uncon.under.txt", sep = ""))
write.table(result.under.mf, paste(outDir, "/go_molecularFunction_uncon.under.txt", sep = ""))
write.table(result.under.bp, paste(outDir, "/go_biologicalProcess_con.under.txt", sep = ""))
write.table(result.under.cc, paste(outDir, "/go_cellularComponent_con.under.txt", sep = ""))
write.table(result.under.mf, paste(outDir, "/go_molecularFunction_con.under.txt", sep = ""))

