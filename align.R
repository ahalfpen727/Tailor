
source("http://bioconductor.org/biocLite.R")


options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number, prefixed by #
## options(echo=TRUE)

# We need to keep track of versions used, especially packages
sessionInfo()

# dir = "~/Documents/UMB/Riley_Lab/HiSeq_data/jose/"
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

library('Rsubread')

readfile1 = paste0(inFile, "_R1.fastq.gz")
readfile2 = paste0(inFile, "_R2.fastq.gz")

print(paste0("readfile1=",readfile1))
print(paste0("readfile2=",readfile2))

## readfile1=c(
##     "createLinks_results_default/Sample_10211_R1.fastq.gz",
##     "createLinks_results_default/Sample_10212_R1.fastq.gz",
##     "createLinks_results_default/Sample_10213_R1.fastq.gz",
##     "createLinks_results_default/Sample_10214_R1.fastq.gz",
##     "createLinks_results_default/Sample_10215_R1.fastq.gz",
##     "createLinks_results_default/Sample_1021B1_R1.fastq.gz",
##     "createLinks_results_default/Sample_1021B2_R1.fastq.gz",
##     "createLinks_results_default/Sample_1021B3_R1.fastq.gz",
##     "createLinks_results_default/Sample_1021B4_R1.fastq.gz",
##     "createLinks_results_default/Sample_1021B5_R1.fastq.gz",
##     "createLinks_results_default/Sample_A1_R1.fastq.gz",
##     "createLinks_results_default/Sample_A2_R1.fastq.gz",
##     "createLinks_results_default/Sample_A3_R1.fastq.gz",
##     "createLinks_results_default/Sample_A4_R1.fastq.gz",
##     "createLinks_results_default/Sample_A5_R1.fastq.gz",
##     "createLinks_results_default/Sample_AB1_R1.fastq.gz",
##     "createLinks_results_default/Sample_AB2_R1.fastq.gz",
##     "createLinks_results_default/Sample_AB3_R1.fastq.gz",
##     "createLinks_results_default/Sample_AB4_R1.fastq.gz",
##     "createLinks_results_default/Sample_AB5_R1.fastq.gz"
## )

## readfile2=c(
##     "createLinks_results_default/Sample_10211_R2.fastq.gz",
##     "createLinks_results_default/Sample_10212_R2.fastq.gz",
##     "createLinks_results_default/Sample_10213_R2.fastq.gz",
##     "createLinks_results_default/Sample_10214_R2.fastq.gz",
##     "createLinks_results_default/Sample_10215_R2.fastq.gz",
##     "createLinks_results_default/Sample_1021B1_R2.fastq.gz",
##     "createLinks_results_default/Sample_1021B2_R2.fastq.gz",
##     "createLinks_results_default/Sample_1021B3_R2.fastq.gz",
##     "createLinks_results_default/Sample_1021B4_R2.fastq.gz",
##     "createLinks_results_default/Sample_1021B5_R2.fastq.gz",
##     "createLinks_results_default/Sample_A1_R2.fastq.gz",
##     "createLinks_results_default/Sample_A2_R2.fastq.gz",
##     "createLinks_results_default/Sample_A3_R2.fastq.gz",
##     "createLinks_results_default/Sample_A4_R2.fastq.gz",
##     "createLinks_results_default/Sample_A5_R2.fastq.gz",
##     "createLinks_results_default/Sample_AB1_R2.fastq.gz",
##     "createLinks_results_default/Sample_AB2_R2.fastq.gz",
##     "createLinks_results_default/Sample_AB3_R2.fastq.gz",
##     "createLinks_results_default/Sample_AB4_R2.fastq.gz",
##     "createLinks_results_default/Sample_AB5_R2.fastq.gz"
## )

output_format="BAM"

align(
                                        # index for reference sequences
    index="1021_Sm_genome",
                                        # input reads and output
    readfile1=readfile1,

    readfile2=readfile2,

    type="rna",
    input_format="gzFASTQ",
    output_format="BAM",
    output_file=paste(as.character(readfile1),"subread",output_format,sep="."),
    # output_file="Rsubread.alignment.BAM",
                                        # offset value added to Phred quality scores of read bases
    phredOffset=33,
                                        # thresholds for mapping
    nsubreads=10,
    TH1=3,
    TH2=1,
    maxMismatches=3,
                                        # unique mapping and multi-mapping
    unique=TRUE,
    nBestLocations=1,
                                        # indel detection
    indels=5,
    complexIndels=FALSE,
                                        # read trimming
    nTrim5=0,
    nTrim3=0,
                                        # distance and orientation of paired end reads
    minFragLength=50,
    maxFragLength=600,
    PE_orientation="fr",
                                        # number of CPU threads
    nthreads=1,
                                        # read group
    readGroupID=NULL,
    readGroup=NULL,
                                        # color space reads
    color2base=FALSE,
                                        # dynamic programming
    DP_GapOpenPenalty=-1,
    DP_GapExtPenalty=0,
    DP_MismatchPenalty=0,
    DP_MatchScore=2,
                                        # detect structural variants
    detectSV=FALSE
                                        # gene annotation
    ## useAnnotation=TRUE,
    ## annot.inbuilt=NULL,
    ## annot.ext="Smeliloti_1021_refseq_bowtie_assembly/1021_genome.gtf",
    ## isGTF=TRUE,
    ## GTF.featureType="CDS",
    ## GTF.attrType="Parent",
    ## chrAliases=NULL
)
