#!/bin/bash

#----NOTE-------------------------------------------------------------------------------------------------------------------------------------------
# Bash does NOT accept spaces between the equal sign and the variable name NOR between the equal sign and the assigned value!
#
#    Example 1: variable=value     <-- Yes!
#    Example 2: variable = value   <-- No!
#    Example 3: variable= value    <-- No!
#    Example 4: variable =value    <-- No!
#---------------------------------------------------------------------------------------------------------------------------------------------------


#----PATH PARAMETERS--------------------------------------------------------------------------------------------------------------------------------
TAILOR=/project/umb_triley/cpct/rna-seq/tailor
RSYNC=rsync_results
#RSYNC=/project/umb_cpct/data/jose/hiseq/12.01.2014/141201_D00345_0037_BHBATBADXX/Data/Intensities/BaseCalls
BCL2FASTQ=bcl2fastq_results
MERGE_FASTQ_FILES=mergeFastqFiles_results
LINKS=createLinks_results
PRE_TRIM_FASTQC=preTrim_fastqc_results
TRIM=trim_results
POST_TRIM_FASTQC=postTrim_fastqc_results
TOPHAT=tophat_results
CUFFLINKS=cufflinks_results
CUFFMERGE=cuffmerge_results
CUFFCOMPARE=cuffcompare_results
CUFFCOMPARE_MERGE=cuffcompare_merge_results
CUFFQUANT=cuffquant_results
CUFFNORM=cuffnorm_results
CUFFDIFF=cuffdiff_results
LEAVEOUT=leaveout_results
GENENAMES=diffGeneNames_results
CUMMERBUND=cummeRbund_results
SAMTOOLS=samtools_results
BAMTOOLS=bamtools_results
BIONET=bionet_results
BOWTIE=bowtie_results
PATHVIEW=pathview_results
GO=go_results
MATRIX=matrix_results
CLASSIFICATION=classification_results

#--GENOME PATHS---------------------------------------------------------------------------
# Using Tailor requires a genome fasta file at the very least. The variables below tell the tools where to find the genome fasta file.  From the genome fasta file a bowtie-index and reference annotation file can be extracted, enabling the tophat-mapping step as well as all the steps that follow tophat. 

GRCH38_BOWTIE2_INDEX=/project/umb_triley/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome
GRCH38_GTF=/project/umb_triley/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
GRCH38_FASTADIR=/project/umb_triley/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes
GRCH38_MULTIFASTA=/project/umb_triley/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
GRCH38_FAIDX=/project/umb_triley/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa.fai

HG19_BOWTIE2_INDEX=/project/umb_triley/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
HG19_GTF=/project/umb_triley/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
HG19_MULTIFASTA=/project/umb_triley/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
HG19_FAIDX=/project/umb_triley/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai
HG19_FASTADIR=/project/umb_triley/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
#HG19_RIBO_BOWTIE_INDEX=/project/umb_triley/cpct/rna-seq/tailor/indexFiles/hg19_ribosomal_index                                                                                     #HG19_TRANSCRIPTOME=/project/umb_triley/cpct/rna-seq/urine1/hg38_tophat_transcriptome     #HG19_MULTIFASTA=/project/umb_triley/cpct/rna-seq/urine1/hg19/tophat_transcriptome_HG19/hg19_perfect.fa                                                                             #HG19_MASK=/project/umb_triley/cpct/rna-seq/drew_mods/hg19/UCSC_HG19/Annotations/phiX.fa  # Other fasta files (RefSeq and Ensembl) in /project/umb_triley/genomes/Homo_sapiens/UCSC # /project/umb_triley/genomes/Homo_sapiens/UCSC/ensembl.ucsc.1000up.1000down.fasta        # /project/umb_triley/genomes/Homo_sapiens/UCSC/refseq.ucsc.1000up.1000down.fasta         


# Viral Genome Spike in for Illumina calibration and may be useful for normalizing/quantifying gene expression
#PHIX_BOWTIE_INDEX=/project/umb_triley/cpct/rna-seq/tailor/indexFiles/phiX174_bowtie_index
#PHIX_MULTIFASTA=/project/umb_triley/cpct/rna-seq/tailor/indexFiles/phiX174
#PHIX_FASTADIR=/project/umb_triley/cpct/rna-seq/tailor/indexFiles/
#PHIX_GTF=/project/umb_triley/cpct/rna-seq/tailor/indexFiles/

#--------------------------------------------------------------------------------------------------------------------------------------------

FASTADIR=$HG19_FASTADIR
MULTIFASTA=$HG19_MULTIFASTA
BOWTIE_INDEX=$HG19_BOWTIE2_INDEX
REFERENCE_GTF=$HG19_GTF
FASTA_INDEX=$HG19_FAIDX

#FASTADIR=$GRCH38_FASTADIR
#MULTIFASTA=$GRCH38_MULTIFASTA
#BOWTIE_INDEX=$GRCH38_BOWTIE2_INDEX
#REFERENCE_GTF=$GRCH38_GTF
#FASTA_INDEX=$GRCH38_FAIDX

JOBS_OUT_DIR=jobs.out
JOBS_SCRIPTS_DIR=jobs.scripts
#--------------------------------------------------------------------------------------------------------------------------------------------

#----LABEL PARAMETER THAT DIFFERENTIATES THIS RUN OF THE PIPELINE FROM OTHER RUNS-------------------------------------------------------------
#----IF THE LABEL IS NOT "_DEFAULT" THEN YOU NEED TO COMMENT ABOUT WHICH PARAMETERS WERE CHANGED----------------------------------------------

DEFAULT_LABEL="_default"
GENOME_LABEL="_hg19"
#GENOME_LABEL="_grch38"

# for differentiating tophat runs with experimental parameters it is advised that an additional label be used to prevent over writes

# GRCH38
#INPUT_LABEL="_hg19_default"
#OUTPUT_LABEL="_hg19_default"
INPUT_LABEL="_grch38_default"
OUTPUT_LABEL="_grch38_default_no_N_no_u"
#OUTPUT_LABEL="_grch38_gtf_only"

# for differentiating cufflinks runs
# OUTPUT_LABEL=".$(date +%F_%R)"  # a date-time stamp, example: ".2015-10-16_19:25"
#---------------------------------------------------------------------------------------------------------------------------------------------------

#----JOB SUBMISSION PARAMETERS---------------------------------------------------------------------
SCHEDULER=lsf
EXTRA_SUBMIT_PARAMETERS="-R (hname!='ghpcc-sgi')"
#--------------------------------------------------------------------------------------------------

#----XVFB-RUN PARAMETERS---------------------------------------------------------------------------
XVFB_RUN="xvfb-run -a -n 1 -s \"-screen 0 1600x1200x24\""
#--------------------------------------------------------------------------------------------------


#----RSYNC PARAMETERS--------------------------------------------------------------------------------------------------------------------------
RSYNC_OUTPUT=${RSYNC}${DEFAULT_LABEL}
#SOURCE_TYPE="*.bcl"
#SOURCE_DIR="/full/path/on/source/hostname"
#SOURCE_HOST="papabear.umb.edu"
#SOURCE_USER="todd"
#SOURCE_PASS_FILE="~/${SOURCE_HOST}.pass"  # create this file that contains only the password and then "chmod 600" it SO THAT NO ONE ELSE CAN VIEW IT!
# WE NEVER PASS THE PASSWORD ON THE COMMAND LINE OR AS AN ENVIRONMENT VARIABLE! THEREFORE, PEOPLE CANNOT SEE IT!!
# EXTRA_RSYNC_PARAMETERS="-auvz"


#----GENOME RSYNC PARAMETERS--------------------------------------------------------------------------------------------------------------------------                     
EXTRA_RSYNC_PARAMETERS=""  #"-auvZ"
#SOURCE_TYPE="*"
#SOURCE_DIR="/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens_Ensembl_GRCh38.tar.gz" 
# ftp://ftp.ensembl.org/../pub/current_fasta/homo_sapiens/dna/
#SOURCE_HOST="ussd-ftp.illumina.com"
#SOURCE_USER="igenome"   # G3nom3s4u
#SOURCE_PASS_FILE="~/${SOURCE_HOST}.pass"  # This is a public key generated by illumina
#location: ftp.broadinstitute.org/bundle
#username: gsapubftp-anonymous
#password:
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_embl/homo_sapiens .
# ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
#----------------------------------------------------------------------------------------------------------------------------------------------

#----BCL2FASTQ PARAMETERS---------------------------------------------------------------------------------------------------------------------
BCL2FASTQ_INPUT=${RSYNC}${DEFAULT_LABEL}
BCL2FASTQ_OUTPUT=$BCL2FASTQ${DEFAULT_LABEL}
SAMPLE_SHEETS=("SampleSheet.csv")
# SAMPLE_SHEETS=("SampleSheet.filledIndex.csv")
# SAMPLE_SHEETS=("SampleSheet.filledIndex.lanes1-4.csv")
# SAMPLE_SHEETS=("SampleSheet.filledIndex.lanes1-4.csv" "SampleSheet.filledIndex.lanes5-8.csv")
BASES_MASK="y51,I6,y51"
EXTRA_BCL2FASTQ_PARAMETERS="--fastq-cluster-count 0 --no-eamss --force  --no-lane-splitting"
#----------------------------------------------------------------------------------------------------------------------------------------------

#----MERGE FASTQ FILES PARAMETERS-------------------------------------------------------------------------------------------------------------
MERGE_INPUT=${BCL2FASTQ}${DEFAULT_LABEL}
MERGE_OUTPUT=${MERGE_FASTQ_FILES}${DEFAULT_LABEL}
EXTRA_MERGE_PARAMETERS=""
MERGE_FIND_INPUT_PATTERN='*'
MERGE_REGEX='(\S+)_(\S+)_(\S+)_(\S+)_(\S+)\.fastq\.gz'  # Group0=everything, Group1=sample, Group2=barCode, Group3=lane
SAMPLE_INDEX=1
BARCODE_INDEX=2
LANE_INDEX=3
R1_GLOB='*${barcode}*_R1_*.fastq.gz'
R1_DESTINATION='${sample}_${barcode}_R1.fastq.gz'
R2_GLOB='*${barcode}*_R2_*.fastq.gz'
R2_DESTINATION='${sample}_${barcode}_R2.fastq.gz'
#-----------------------------------------------------------------------------------------------------------------------------------------------

#----LINK PARAMETERS---------------------------------------------------------------------------------------------------------------------------

#  The links should have VERY SHORT meaningful names constructed from the following keywords:

#  CTRL = control, CASE = case, EXP  = experiment, SAMP = sample,  SYMP = symptomatic,  LN   = lane, MPLX = multiplex

#-------------------------------------------------------------------------------------------------------------------------------------------- 
LINKS_INPUT=${MERGE_FASTQ_FILES}${DEFAULT_LABEL}
LINKS_OUTPUT=${LINKS}${DEFAULT_LABEL}
LINKS_FILE=createLinks.table  # tab separated table where 1st column is the fastq source file and the 2nd column is the destination link
EXTRA_LINK_PARAMETERS=""
#----------------------------------------------------------------------------------------------------------------------------------------------

#----FASTQC PARAMETERS------------------------------------------------------------------------------------------------------------------------------
FASTQC_INPUT=${LINKS}${DEFAULT_LABEL}
FASTQC_OUTPUT=${PRE_TRIM_FASTQC}${INPUT_LABEL0}
FASTQC_FIND_INPUT_PATTERN='*.fastq.gz'
FASTQC_INPUT_REGEX="(.+)_R(.+)\.fastq\.gz"  # Group1=sample, Group2=pairedNumber
EXTRA_FASTQC_PARAMETERS=""
#--------------------------------------------------------------------------------------------------------------------------------------------

#----TRIM PARAMETERS--------------------------------------------------------------------------------------------------------------------------
TRIM_INPUT=${LINKS}${DEFAULT_LABEL}
TRIM_OUTPUT_FASTQC="${POST_TRIM_FASTQC}${DEFAULT_LABEL}"
TRIM_OUTPUT_TRIM="${TRIM}${DEFAULT_LABEL}"
TRIM_FIND_INPUT_PATTERN='*.fastq.gz'
TRIM_INPUT_REGEX="(.+)_R(.+)\.fastq\.gz"  # Group1=sample, Group2=pairedNumber
# EXTRA_TRIM_PARAMETERS='--clip_R1 8 --paired --retain_unpaired ${R1file} ${R2file}'
EXTRA_TRIM_PARAMETERS='--paired --retain_unpaired ${R1file} ${R2file}'
#-----------------------------------------------------------------------------

#----TOPHAT PARAMETERS------------------------------------------------------------------------------------------------------------------------------
TOPHAT_INPUT=${TRIM}${DEFAULT_LABEL}
TOPHAT_OUTPUT=${TOPHAT}${GENOME_LABEL}${DEFAULT_LABEL}
TOPHAT_FIND_INPUT_PATTERN='*val*.fq.gz'
TOPHAT_INPUT_REGEX="(.+)_R1"
# GTF_GUIDED analyses perform gtf assisted de novo discovery. 
EXTRA_TOPHAT_PARAMETERS='-g 1 -G ${REFERENCE_GTF} --coverage-search --microexon-search --fusion-search --library-type fr-firststrand'

#EXTRA_TOPHAT_PARAMETERS='-g 1 -G ${REFERENCE_GTF} --coverage-search --microexon-search --fusion-search' 
# EXTRA_TOPHAT_PARAMETERS='-g 1 -G ${CUFFCOMPARE_MERGE}${OUTPUT_LABEL}/cuffcmp.combined.gtf --coverage-search --microexon-search --fusion-search --library-type fr-firststrand'
# DEFAULT and GTF_ONLY analyses increase concordant mapping + exclude unannotated splice junc (not in the gtf gff file). 
# '-g 1  -G ${REFERENCE_GTF} --no-novel-juncs'

# to determine if tophat (or any step) has completed and which samples failed use:
# > grep -rc "Successfully completed" ./jobs.out/*out
# to determine mapping rates use:
# > grep -r concordant ./Sample_*/align_summary.txt >> concordant_mapping                                         
# > grep -r overall ./Sample_*/align_summary.txt >> overall_mapping                                                        
# > grep -r multiple ./Sample_*/align_summart.txt >> multiple_mapping 
#----- DESCRIPTION OF SOME OF TOPHAT's OPTIONS----------------
# BE AWARE of additional parameters that are hard coded into Tophat and Cuffdiff steps (listed at beginning of the settings file). 
# The following options can be included in the ${EXTRA_TOPHAT_PARAMETERS} variable
# "--coverage-search 
# "--microexon-search --fusion-search# searches for reads incident upon microexons and gene fusions.           
# "--transcriptome-only"  # Only aligns the reads to the transciptome and report only those mappings as genomic mappings
# "--read-realign-edit-dist=0" # best alignment after realignment (time intensive)                                       
# " -g 10" # Maximum number of multi-hits (multiple allowable mappings) to 10                                            
#EXTRA_TOPHAT_PARAMETERS='-G ${REFERENCE_GTF} --transcriptome-index ${ANNOTATION_INDEX}' # builds a transcriptome (exome) from the wholegenome.fa and GTF file       
#------------------------------------------------------------------------------


#---SAMTOOLS PARAMETERS-----------------------------------------------------------------------------------------------------------------------------
SAMTOOLS_INPUT=${TOPHAT}${OUTPUT_LABEL}
#SAMTOOLS_OUTPUT="${SAMTOOLS}${OUTPUT_LABEL}" # add samtools label _toolused to differentiate and create outdir                        
SAMTOOLS_OUTPUT="${SAMTOOLS}${OUTPUT_LABEL}_flagstat"
# produces mapping statistics similar to unix command search but much lower limit for acceptable intronic length 
#SAMTOOLS_OUTPUT="${SAMTOOLS}${OUTPUT_LABEL}_depth"                                                                                  
#SAMTOOLS_OUTPUT="${SAMTOOLS}${OUTPUT_LABEL}_mpileup"  
# depth and mpileup perform coverage assessment and SNP analysis of the input SAM file
# REGION_BED=/project/umb_triley/cpct/rna-seq/tailor/indexFiles/chr10_region.bed
# REGION=chr10:44865600-44880545                            
# REGION-> #Genomic region chr##:####-####
SAMTOOLS_FIND_INPUT_PATTERN='*_out'
SAMTOOLS_INPUT_REGEX="(.+)_out"   # anything that ends with "_out"
EXTRA_SAMTOOLS_PARAMETERS=""
#-----------------------------------------------------------------------------------------------------------------------------------------------

#----CUFFQUANT_PARAMETERS---------------------------------------------------------------------------------------------------------------------------
CUFFQUANT_INPUT="${TOPHAT}${OUTPUT_LABEL}"
CUFFQUANT_OUTPUT="${CUFFQUANT}${OUTPUT_LABEL}"
CUFFQUANT_INPUT_REGEX="(.+)_out"   # anything that ends with "_out"
CUFFQUANT_FIND_INPUT_PATTERN='*_out'
CUFFUANT_FIND_INPUT_FILE='input.list.txt'
CUFFQUANT_ALIGNMENT_FILE="accepted_hits.bam"
EXTRA_CUFFQUANT_PARAMETERS="-v ${REFERENCE_GTF} -b ${MULTIFASTA}"
#EXTRA_CUFFQUANT_PARAMETERS="-v ${CUFFCOMPARE_MERGE}${OUTPUT_LABEL}/cuffcmp.combined.gtf -b ${MULTIFASTA}"
# EXTRA_CUFFQUANT_PARAMETERS="-v ${CUFFMERGE}${OUTPUT_LABEL1}/merged.gtf -b ${MULTIFASTA}"
# Cuffquant analyzes and outputs expression values for transcripts in the bam files against a reference gtf file (cuffcmp.combined.gtf )
# may use either cufflinks --> cuffmerge --> cuffcompare cuffquant/cuffdiff    OR  cufflinks --> cuffcomparemerge ---> cuffquant/cuffdiff
# OR cufflinks --> cuffmerge --> cuffquant/cuffdiff
#-----------------------------------------------------------------------------------------------------------------------------------------------

#----CUFFDIFF_PARAMETERS------------------------------------------------------------------------------------------------------------------------
CUFFDIFF_INPUT_GTF="${REFERENCE_GTF}"  # simpler (no novel genes included). Requires reference genome and reference annotation in GFF format

# cufflinks paper protocol parameters (fully de novo)
#CUFFDIFF_INPUT_GTF="${CUFFMERGE}${OUTPUT_LABEL}/merged.gtf"  

# reference annnotation guided parameter
#CUFFDIFF_INPUT_GTF="${CUFFCOMPARE_MERGE}${OUTPUT_LABEL}/cuffcmp.combined.gtf"

#CUFFDIFF_INPUT_ALIGNMENTS="${TOPHAT}${INPUT_LABEL}"  # This option can be used during default (or any analysis) but it is much slower and may fail
#CUFFDIFF_ALIGNMENT_FILE="accepted_hits.bam"

CUFFDIFF_INPUT_ALIGNMENTS="${CUFFQUANT}${OUTPUT_LABEL}"
CUFFDIFF_ALIGNMENT_FILE="abundances.cxb"
#  When cuffquant is included, the cuffquant output, cxb files, are used as the input alignments.
#  This modification reducues compuational load and CPU time for the cuffdiff step.

CUFFDIFF_OUTPUT="${CUFFDIFF}${OUTPUT_LABEL}"
CUFFDIFF_FOLDCHANGE_OUTPUT='${experiment2Name}-over-${experiment1Name}'
EXTRA_CUFFDIFF_PARAMETERS='-L ${experiment1Name},${experiment2Name} -b ${MULTIFASTA}' 
#EXTRA_CUFFDIFF_PARAMETERS='-N -u -L ${experiment1Name},${experiment2Name} -b ${MULTIFASTA}' # -M $HG19_RIBO_GTF'       

# Group 0: CTRL
experimentNames[0]="CTRL";
declare -a experimentGroups0=( \
    Sample_3 \
    Sample_5 \
    Sample_7 \
    Sample_9 \
    Sample_11 \
    Sample_13 \
    Sample_15 \
    Sample_17 \
    Sample_19 \
    );

# Group 1: LUTS
experimentNames[1]="LUTS";
declare -a experimentGroups1=( \
    Sample_4 \
    Sample_6 \
    Sample_8 \
    Sample_10 \
    Sample_12 \
    Sample_14 \
    Sample_16 \
    Sample_18 \
    Sample_20 \
    );

# The following varibale set the comparisons between groups or conditions
NUM_COMPARISONS='1'
# Comparison 1: LUTS over CTRL
 declare -a cuffdiff0=(0 1)  # cuffdiff(group0 group1) calculates group1-over-group0 expression fold-change

# Comparison 2: GROUP2 over GROUP0
# declare -a cuffdiff1=(0 2)  # cuffdiff(group0 group2) calculates group2-over-group0 expression fold-change

# Comparison 3: GROUP3 over GROUP0
# declare -a cuffdiff2=(0 3)  # cuffdiff(group0 group3) calculates group3-over-group0 expression fold-change
#--------------------------------------------------------------------------------------------------------------------------------

#----CUFFNORM_PARAMETERS--------------------------------------------------------------------------------------------------------------------
#
# Normalizes cufflinks fpkm values to counts in a table format or in cuffdiff format to be used by count based tools like edgeR and DEseq
#
#----------------------------------------------------------------------------------------------------------                                      
CUFFNORM_INPUT="${CUFFQUANT}${OUTPUT_LABEL}"
CUFFNORM_OUTPUT="${CUFFNORM}${OUTPUT_LABEL}"

# default, 3-step, simple analysis of known and annotated genes (no de novo discovery)
CUFFNORM_INPUT_GTF="${REFERENCE_GTF}"  

# Cufflinks Paper Protocol (Fully de novo)
#CUFFNORM_INPUT_GTF="${CUFFMERGE}${OUTPUT_LABEL}/merged.gtf"  

# Gtf guided/ Reference Annotation Guided Parameters
#CUFFNORM_INPUT_GTF="${CUFFCOMPARE_MERGE}${OUTPUT_LABEL}/cuffcmp.combined.gtf"

#CUFFNORM_INPUT_ALIGNMENTS="${TOPHAT}${INPUT_LABEL}"  # This option can be used during default (or any analysis) but it is much slower and may fail
#CUFFNORM_ALIGNMENT_FILE="accepted_hits.bam"

CUFFNORM_INPUT_ALIGNMENTS="${CUFFQUANT}${OUTPUT_LABEL}"
CUFFNORM_ALIGNMENT_FILE="abundances.cxb"
#  When cuffquant is included, the cuffquant output, cxb files, are used as the input alignments.
#  This modification reducues compuational load and CPU time for the cuffnorm step.

CUFFNORM_FOLDCHANGE_OUTPUT='${experiment2Name}-over-${experiment1Name}'

# The following varibale set the comparisons between groups or conditions
NUM_COMPARISONS='1'
# Comparison 1: GROUP1 over CTRL
 declare -a cuffnorm0=(0 1)  # cuffnorm(group0 group1) calculates group1-over-group0 expression fold-change

# Comparison 2: GROUP2 over GROUP0
# declare -a cuffnorm1=(0 2)  # cuffnorm(group0 group2) calculates group2-over-group0 expression fold-change

# Comparison 3: GROUP3 over GROUP0
# declare -a cuffnorm2=(0 3)  # cuffnorm(group0 group3) calculates group3-over-group0 expression fold-change

EXTRA_CUFFNORM_PARAMETERS='${CUFFNORM_INPUT_GTF} -output-format cuffdiff -L ${experiment1Name},${experiment2Name}'

#----------------------------------------------------------------------------------------------------------------------------------------

#-----------LAUNCHR_PARAMETERS-------------------------------------------------------------------------------------------------------------------------
# launches the cummerbund, bionet, go, and pathview R scripts via the launchr shell scripts
# CummeRbund.R uses the launchR.eachDiffDir.sh script, while GO.R, Pathview.R and Bionet.R use the launchR.eachDiffFile.sh script
#----------------------------------------------------------------------------------------------------------------------------------------
# To run the R analysis steps on CUFFNORM data, 
# change ${LAUNCHR_EACHDIFFFILE_FIND_INPUT_PATTERN}='*fpkm_tracking'
# If you are running cummeRbund on a de novo cuffdiff output, you may need to increase the expected memory usage with the variable below
# The GTF_GUIDED and DEFAULT analyses should complete with 32 GB, the DE_NOVO analysis may need upward of 70 GB
# MEMORY=80192
#----------------------------------------------------------------------------------------------------------------------------------------

#LAUNCHR_INPUT="${CUFFNORM}${OUTPUT_LABEL}"
LAUNCHR_INPUT="${CUFFDIFF}${OUTPUT_LABEL}"
LAUNCHR_OUTPUT='${1}_results${OUTPUT_LABEL}'
# The INPUT PATTERN below is modified to take cuffdiff output as input 
LAUNCHR_EACHDIFFFILE_FIND_INPUT_PATTERN='*.diff'
# These LaunchR variables are for the CummeRbund step
LAUNCHR_EACHDIFFDIR_FIND_INPUT_DIR='${INPUT}'
LAUNCHR_EACHDIFFDIR_FIND_INPUT_PATTERN='*-over-*'
# These LaunchR variables are for the Pathview, GO, and Bionet (needs work)
LAUNCHR_EACHDIFFFILE_FIND_INPUT_DIR='${INPUT}/*-over-*'
# The variable below has been modified to take cuffnorm output as input. Be sure to change the INPUT PATTERN when using cuffdiff as input 
# LAUNCHR_EACHDIFFFILE_FIND_INPUT_PATTERN='*fpkm_tracking'  
LAUNCHR_INPUT_REGEX="(.+)-over-(.+)"


#----------------------------------------------------------------------------------------------------------------------------------------------


#----CUMMERBUND PARAMETERS--------------------------------------------------------------------------------------------------------------------------
# Default Parameters
#--------------------------------------------------------------------------------------------------------------------------------------------

CUMMERBUND_OUTPUT="${CUMMERBUND}${OUTPUT_LABEL}"
CUMMERBUND_INPUT_GTF="${REFERENCE_GTF}"                                                                                               
GENE_LIST="CXCL12,CXCR4,TGFB,MAPK,IL8,CCL5"
RPLOTS="Rplots.pdf"
FPKM_MATRIX="fpkm.matrix.csv"
DIFF_TABLE="DiffTable.csv"
R_GENOME="hg38"

export DIFF_TABLE
export RPLOTS
export FPKM_MATRIX
export R_GENOME
export CUMMERBUND_INPUT_GTF
export GENE_LIST

# Gtf guided/ Reference Annotation Guided Parameters
#CUMMERBUND_INPUT_GTF="${CUFFCOMPARE_MERGE}${OUTPUT_LABEL}/cuffcmp.combined.gtf"
# Cufflinks paper (fully de novo)
#CUMMERBUND_INPUT_GTF="${CUFFMERGE}${OUTPUT_LABEL}/merged.gtf"

#---------------------------------------------------------------------------------------------------------------------------------------------------

#----BIONET PARAMETERS------------------------------------------------------------------------------------------------------------------------------
 BIONET_OUTPUT=${BIONET}${OUTPUT_LABEL}
 BIONET_GENE_ALPHA_VALUE=0.5
 BIONET_SUBNET_ALPHA_VALUE=0.1
 BIONET_NUMBER_NETWORKS=20

 export BIONET_GENE_ALPHA_VALUE
 export BIONET_SUBNET_ALPHA_VALUE
 export BIONET_NUMBER_NETWORKS

#----GO PARAMETERS----------------------------------------------------------------------------------------------------------------------------------
GO_OUTPUT=${GO}${OUTPUT_LABEL}
GO_DIFF_EXPRESSED_ALPHA_VALUE=0.1   # Alpha value for classifying genes as differentially expressed
GO_HYPER_GEO_ALPHA_VALUE=0.01   #Alpha value for enrichment of differentially expressed genes is a category

export GO_DIFF_EXPRESSED_ALPHA_VALUE
export GO_HYPER_GEO_ALPHA_VALUE
#---------------------------------------------------------------------------------------------------------------------------------------------------


#----PATHVIEW PARAMETERS----------------------------------------------------------------------------------------------------------------------------
PATHVIEW_OUTPUT=${PATHVIEW}${OUTPUT_LABEL}
PATHVIEW_SPECIES="human"   # pathview common.name ("human" or "mouse", etc.)
PATHVIEW_ALPHA_DOWN=0.2  # Alpha value for genes that are under expressed in the case
PATHVIEW_ALPHA_UP=0.2   # Alpha value for genes that are over expressed in the case
PATHVIEW_NUMBER_PATHWAYS_DOWN=30
PATHVIEW_NUMBER_PATHWAYS_UP=30
PATHVIEW_PATHWAYS_LIST=""

export PATHVIEW_SPECIES
export PATHVIEW_ALPHA_DOWN
export PATHVIEW_ALPHA_UP
export PATHVIEW_NUMBER_PATHWAYS_DOWN
export PATHVIEW_NUMBER_PATHWAYS_UP
export PATHVIEW_PATHWAYS_LIST
#-----------------------------------------------------------------------------------------------------------------------------------------------
#----MATRIX PARAMETERS-------------------------------------------------------------------------------------------------------------
MATRIX_OUTPUT=${MATRIX}${OUTPUT_LABEL}
MATRIX_INPUT=${CUFFDIFF}${OUTPUT_LABEL}/
MATRIX_VALUE_X=""
MATRIX_VALUE_n=""
 
#-----------------------------------------------------------------------------------------------------------------------------------------------


#----CLASSIFICATION PARAMETERS----------------------------------------------------------------------------------------------------------------------------


### classification input directory which is equal to output of the Cummerbund
LAUNCHR_CLASSIFICATION_INPUT="${CUMMERBUND_OUTPUT}/*-over-*/${FPKM_MATRIX}"


CLASSIFICATION_OUTPUT=${CLASSIFICATION}${OUTPUT_LABEL}
METHOD=(lasso elasticNet ridge)
FOLD=10                # -1 = leave-one-out
NESTED="NO"              # YES/NO
Iterations=100
BootstrapIterations=0  # 0 = no boostrap
TRANSFORM="NA"           # (log, log2. log10, exp)
MAXBIOMARKERS=500
SCALE="none" #row, col
 
export METHOD
export FOLD
export NESTED
export Iterations
export BootstrapIterations
export TRANSFORM
export MAXBIOMARKERS
export SCALE

ALL_RESULTS_DIRS=( \
    ${RSYNC_OUTPUT} \
    ${BCL2FASTQ_OUTPUT} \
    ${MERGE_FASTQ_FILES_OUTPUT} \
    ${LINKS_OUTPUT} \
    ${FASTQC_OUTPUT} \
    ${TRIM_OUTPUT} \
    ${TOPHAT_OUTPUT} \
    ${CUFFLINKS_OUTPUT} \
    ${CUFFMERGE_OUTPUT} \
    ${CUFFCOMPARE_OUTPUT} \
    ${CUFFCOMPARGE_MERGE_OUTPUT} \
    ${CUFFQUANT_OUTPUT} \
    ${CUFFNORM_OUTPUT} \
    ${SAMTOOLS_OUTPUT} \
    ${BAMTOOLS_OUTPUT} \
    ${CUFFDIFF_OUTPUT} \
    ${LEAVEOUT_OUTPUT} \
    ${GENENAMES_OUTPUT} \
    ${CUMMERBUND_OUTPUT} \
    ${BIONET_OUTPUT} \
    ${BOWTIE2INSPECT_OUTPUT} \
    ${BOWTIE2BUILD_OUTPUT} \
    ${PATHVIEW_OUTPUT} \
    ${GO_OUTPUT} \
    ${MATRIX_OUTPUT} \
    ${CLASSIFICATION_OUTPUT}\
    );
