## source('openGraphSaveGraph.r');
options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number, prefixed by #

### R code from vignette source 'vignettes/cummeRbund/inst/doc/cummeRbund-manual.Rnw'
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

###################################################################################################
## To reverse the direction of a particular column, the method depends on the data type:

##     Numbers: put a - in front of the variable name, e.g. df[ order(-df$weight), ].
##     Factors: convert to integer and put a - in front of the variable name, e.g. df[ order(-xtfrm(df$size)), ].
##     Characters: there isn’t a simple way to do this. One method is to convert to a factor first and then sort as above.
###################################################################################################

diffFiles = list.files(path=inDir, pattern="\\.diff$")

options(java.parameters = "-Xmx16000m")

library(xlsx)
library(tools)

wb = createWorkbook()

# loop over all *.diff files
for (i in 1:length(diffFiles)) {
    label = file_path_sans_ext(diffFiles[i])
    label

    # reaads in the dataframe "data"
    data = read.table(paste(inDir, diffFiles[i], sep="/"), sep="\t", header = TRUE)

    if (grepl("_exp", label)) { # does diffFile ends with "_exp"?

        # add fold-change column
        data$fold_change = 2 ** data[, "log2.fold_change."] # as a header, "(" is converted to "."
        # data$fold_change = data[, "log2.fold_change."] # as a header, "(" is converted to "."

        # re-order columns
        #data = data[ , -1]

	#re-order columns
	data = data[ , c(1:9, 15, 10:14)]

        # sort by (1) descending-significant?, then (2) descending-fold-change
        data = data[order(-xtfrm(data$significant), -data$fold_change), ]

    }
    else { # does NOT contain "_exp"

        # sort by (1) descending-significant?, then (3) ascending-q-value
        data = data[order(-xtfrm(data$significant), data$q_value), ]

    }

    sheet1 = createSheet(wb, sheetName = label)
    addDataFrame(data, sheet1)

    #autoSizeColumn(sheet1, 1:14)
    #col.names <- as.character(sheet1)

    # set widths to 20
    #setColumnWidth(sheet1, colIndex=1:ncol(data), colWidth=20)
    
    # autosize column widths
    autoSizeColumn(sheet1, colIndex=1:ncol(data))
    #createFreezePane(sheet, rowSplit, colSplit, startRow=NULL, startColumn=NULL)
    createFreezePane(sheet1, 2, 2, 2, 2)
}

## label = file_path_sans_ext("gene_exp.diff")
## data = read.table(paste(inDir, "gene_exp.diff", sep="/"))
## sheet1 = createSheet(wb, sheetName = label)
## addDataFrame(data, sheet1)


###################################################
### set working directory to output directory
###################################################
setwd(outDir)

#save workbook
saveWorkbook(wb, paste(over, "-over-", under, ".xlsx", sep = ""))

###################################################
### code chunk number 92: session
###################################################
sessionInfo()


