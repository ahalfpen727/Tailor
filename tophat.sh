#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# TOPHAT - aligns reads to the genome and calculates alignments, junctions, insertions, and deletions
#
# This produces the $TOPHAT folder with the following files:
# 	accepted_hits.bam
# 	junctions.bed
# 	insertions.bed
# 	deletions.bed
#
#  - can send you emails when each job starts to execute and when each is finished so
#    that you know when to submit the jobs for the next step
#
#--------------------------------------------------------------------------------------------------------------

#----SETTINGS---------------
if [ -z ${TAILOR_CONFIG+x} ]; then
    echo -e "\nEnvironment variable TAILOR_CONFIG is not set!\n"
    exit
elif [ ! -f "${TAILOR_CONFIG}" ]; then
    echo -e "\nSettings file ${TAILOR_CONFIG} does not exist!\n"
    exit
else
    echo -e "\nUsing settings file ${TAILOR_CONFIG}.\n"
fi

source "${TAILOR_CONFIG}"
#---------------------------

#----JOB SUBMISSION PARAMETERS---------------------------------------------------------------------
#PROCESSORS=6
PROCESSORS=9
MEMORY="7592"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
#MEMORY="5592"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="3:59"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
#DURATION="96:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"        # short = max 4 hours;   long = max 30 days
#QUEUE="long"        # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${TOPHAT_INPUT}"
OUTPUT="${TOPHAT_OUTPUT}"
SCRIPTS=${OUTPUT}/${JOBS_SCRIPTS_DIR}
JOBS=${OUTPUT}/${JOBS_OUT_DIR}
#--------------------------------------------------------------------------------------------------

#----OUTPUT------------------
if [ ! -d ${OUTPUT} ]; then
    mkdir ${OUTPUT}
fi
if [ ! -d ${SCRIPTS} ]; then
    mkdir ${SCRIPTS}
fi
if [ ! -d ${JOBS} ]; then
    mkdir ${JOBS}
fi
if [ ! -d ${ANNOTATION_INDEX} ]; then
    mkdir ${ANNOTATION_INDEX}
fi
#---------------------------

# Save this version of the config file
configFileName=$(basename ${TAILOR_CONFIG})
cp "${TAILOR_CONFIG}" "${OUTPUT}/${configFileName}.$(date +%F_%R)" 

#----GLOB INPUT--------------------------------
findResults=${OUTPUT}/input.list.txt
find ${INPUT} -name "${TOPHAT_FIND_INPUT_PATTERN}" > $findResults
#----------------------------------------------

COMMAND=tophat

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a fastq file

    fileName=$(basename $line)

    if [[ $fileName =~ $TOPHAT_INPUT_REGEX ]]; then

	SAMPLE=${BASH_REMATCH[1]}

	scriptString="mktemp -p ${SCRIPTS}  ${COMMAND}.${SAMPLE}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	chmodString="chmod 777 ${tempScript}"
	echo -e `${chmodString}`

   ### When running Tophat, supply a reference genome annotation file (GTF or GFF).  The transcriptome index and bowtie2 index can then be built ##
   ### with the -G option pointing to the GTF or GFF genome annotation and  the --transcriptome-index option pointing to the directory that will ##
   ### house the newly built transcriptome.  This initial tophat run can be perfomed without performing mapping/aligning to reduce CPU time. ######

	EXTRA_PARAMETERS=$(eval echo "$EXTRA_TOPHAT_PARAMETERS")
	echo -e "\n${EXTRA_PARAMETERS}"

# Cannonical Command
COMMAND_LINE="${COMMAND} -p ${PROCESSORS} ${EXTRA_PARAMETERS} -o ${OUTPUT}/${SAMPLE}_out ${GENOME_INDEX} ${INPUT}/${SAMPLE}_R1_val_1.fq.gz ${INPUT}/${SAMPLE}_R2_val_2.fq.gz"

# The following modification was performed to analyze toxoplasmo data with tailor from an intermediate step
#COMMAND_LINE="${COMMAND} -p ${PROCESSORS} ${EXTRA_PARAMETERS} -o ${OUTPUT}/${SAMPLE}_out ${GENOME_INDEX} ${INPUT}/${SAMPLE}_R1.fastq.gz ${INPUT}/${SAMPLE}_R2.fastq.gz"

	echo -e "source loadModules.sh\n\n" > ${tempScript}
	echo "$COMMAND_LINE" >> ${tempScript}

	echo -e "grep -i multiple ./${OUTPUT}/${SAMPLE}_out/align_summary.txt >> ${OUTPUT}/multiple_mapping.txt\n" >> ${tempScript}
	echo -e "grep -i concordant ./${OUTPUT}/${SAMPLE}_out/align_summary.txt >> ${OUTPUT}/concordant_mapping.txt\n" >> ${tempScript}
	echo -e "grep -i overall ./${OUTPUT}/${SAMPLE}_out/align_summary.txt >> ${OUTPUT}/overall_mapping.txt\n" >> ${tempScript}

	if [ $SCHEDULER == "sge" ]; then
            SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${SAMPLE}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
	else
	    SUBMIT_COMMAND="bsub -q $QUEUE -J ${SAMPLE}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}-${SAMPLE}.%J.out -e ${JOBS}/${COMMAND}-${SAMPLE}.%J.error bash ${tempScript}"
	fi

	date=`date`
	echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
	cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

	echo `${SUBMIT_COMMAND}`
	# rm ${tempScript}
    # else
    #    echo "$line does not match reg.exp."

    fi

done < "$findResults"  # Sends the $findResults file as input to the while-loop
