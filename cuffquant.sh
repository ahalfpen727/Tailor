#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# Cuffquant - quantifies gene expression for each sample and produces a binary cxb file to be used by cuffdiff
#
# In NGS RNA-seq samples, quantitative gene expression data is normalized for total gene/transcript length and the number of sequencing reads, and reported as
#   RPKM: Reads Per Kilobase of exon per Million mapped reads. Used for reporting data based on single-end reads
#   FPKM: Fragments Per Kilobase of exon per Million fragments. Used for reporting data based on paired-end fragments
#
#  - Can send you emails when each job starts to execute and when each is finished so
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
PROCESSORS=3
MEMORY="25192"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="92:00"  # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"       # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${CUFFQUANT_INPUT}"
OUTPUT="${CUFFQUANT_OUTPUT}"
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
#---------------------------

# Save this version of the config file
configFileName=$(basename ${TAILOR_CONFIG})
cp "${TAILOR_CONFIG}" "${OUTPUT}/${configFileName}.$(date +%F_%R)" 

#----GLOB INPUT--------------------------------
 findResults=${OUTPUT}/input.list.txt
#findResults=${OUTPUT}/threeFailed.list.txt
 find ${INPUT} -type d -maxdepth 1 -name "${CUFFQUANT_FIND_INPUT_PATTERN}" > $findResults  # search for only directories and no subdirs
#----------------------------------------------

COMMAND=cuffquant

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a directory that contains a accepted_hit.bam file

    if [[ $line =~ $CUFFQUANT_INPUT_REGEX ]]; then

	i=$(basename ${BASH_REMATCH[1]})

	COMMAND_LINE="${COMMAND} -p ${PROCESSORS} ${EXTRA_CUFFQUANT_PARAMETERS} -o ${OUTPUT}/${i}_out ${INPUT}/${i}_out/${CUFFQUANT_ALIGNMENT_FILE}"

	scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${i}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	echo -e "source loadModules.sh\n\n" > ${tempScript}
	echo "$COMMAND_LINE" >> ${tempScript}
	chmodString="chmod 777 ${tempScript}"
        echo -e `${chmodString}` 

	if [ $SCHEDULER == "sge" ]; then
            SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${i}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
	else
	    SUBMIT_COMMAND="bsub -q $QUEUE -J ${i}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}-${i}.%J.out -e ${JOBS}/${COMMAND}-${i}.%J.error bash ${tempScript}"
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
