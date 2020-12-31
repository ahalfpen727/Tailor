##############################################################################################
# Always version so that the scripts won't break if the "current version" changes!
##############################################################################################

#Unload any software lingering in the environment 
# module purge could replace individual module unloading
module unload fastqc
module unload tophat
module unload cufflinks
module unload bowtie2
module unload samtools
module unload bedtools
module unload sratoolkit
module unload bcl2fastq
module unload trim_galore
module unload perl
module unload R
module unload openssl
module unload python
module unload tablemaker
module unload hisat
module unload bcftools
module unload anaconda

# Initiliaze only the necessary modules
module load python/2.7.5
module load fastqc/0.10.1
module load hisat2/2.0.5
module load tophat/2.0.14
module load cufflinks/2.2.1
module load bowtie2/2-2.1.0
module load samtools/0.0.19
module load sratoolkit/2.3.4-2
module load perl/5.10.1
module load bcl2fastq/1.8.4
module load bcftools/2014_06_02
module load cairo/1.12.16
module load trim_galore/0.3.7
module load bedtools/2.26.0
module load R/3.2.2
module load anaconda3/2019.03

# module load R/3.4.0
# module load R/3.3.1
# module load tablemaker/2.1.1
# module load bcl2fastq/2.17.1.14
# module load bowtie2/2.3.2
# module load samtools/1.4.1
# module load python/2.7.9
# module load python/2.7.9_packages/HTSeq/0.6.1
# module unload python3
# module load fastqc/1.0

source .Rprofile.envrc
