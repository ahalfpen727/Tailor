#!/bin/bash

#-------------------------------------------------
# bcl2fastq - convert BCL files to FASTQ files
#
#  - can send you emails when each job starts to execute and when each is finished so
#    that you know when to submit the jobs for the next step
#
#-------------------------------------------------

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
PROCESSORS=4
MEMORY="8192"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="4:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"       # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS------------------
INPUT=$BCL2FASTQ_INPUT
OUTPUT=$BCL2FASTQ_OUTPUT
SCRIPTS=${OUTPUT}/${JOBS_SCRIPTS_DIR}
JOBS=${OUTPUT}/${JOBS_OUT_DIR}
#---------------------------

#----OUTPUT-----------------
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

########################################################################
# DO NOT SPLIT THE FASTQ FILES!
# ONE FASTQ FILE PER BCL FILE SO THAT IT'S EASY TO KEEP TRACK AND MATCH!
########################################################################

COMMAND=bcl2fastq


for SAMPLE_SHEET in "${SAMPLE_SHEETS[@]}"; do

	COMMAND1=configureBclToFastq.pl
	COMMAND1_LINE="${COMMAND1} --input-dir ${INPUT} --output-dir ${OUTPUT} --sample-sheet ${INPUT}/${SAMPLE_SHEETS} --use-bases-mask ${BASES_MASK}  ${EXTRA_BCL2FASTQ_PARAMETERS}"

	echo -e "SAMPLE_SHEET = ${SAMPLE_SHEET}\n"
	echo -e "COMMAND1_LINE = ${COMMAND1_LINE}\n"

	COMMAND2=make
	COMMAND2_LINE="${COMMAND2} -j $PROCESSORS"

	scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	chmodString="chmod 777 ${tempScript}"
	echo -e `${chmodString}`

	echo -e "source loadModules.sh\n" > ${tempScript}
	echo -e "$COMMAND1_LINE\n"          >> ${tempScript}
	echo -e "cd ${OUTPUT}\n"            >> ${tempScript}
	echo -e "$COMMAND2_LINE\n"          >> ${tempScript}

	if [ $SCHEDULER == "sge" ]; then
	    SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
	else
	    SUBMIT_COMMAND="bsub -K -q $QUEUE -J ${COMMAND} -n ${PROCESSORS} -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}.%J.out -e ${JOBS}/${COMMAND}.%J.error bash ${tempScript}"
	fi

	date=`date`
	echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
	cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

	echo `${SUBMIT_COMMAND}`
# rm ${tempScript}

done

