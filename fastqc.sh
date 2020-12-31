#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# fastqc - perform Quality Control on the fastq files
# This creates 2 directories with output in HTML format for visual inspection in a browser, key statistics and tests are in text file
# The zip files generated can be deleted
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
PROCESSORS=4
MEMORY="8192"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="4:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"        # # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-------------------
INPUT=${FASTQC_INPUT}
OUTPUT=${FASTQC_OUTPUT}
SCRIPTS=${OUTPUT}/${JOBS_SCRIPTS_DIR}
JOBS=${OUTPUT}/${JOBS_OUT_DIR}
#----------------------------

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
find ${INPUT} -name "${FASTQC_FIND_INPUT_PATTERN}" > $findResults
#----------------------------------------------

COMMAND=fastqc

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a fastq file

    fileName=$(basename $line)

    if [[ $fileName =~ $FASTQC_INPUT_REGEX ]]; then

	SAMPLE="${BASH_REMATCH[1]}"
	PAIREDNUMBER="${BASH_REMATCH[2]}"

	# echo "SAMPLE=${SAMPLE}; PAIREDNUMBER=${PAIREDNUMBER}"
	COMMAND_LINE="$COMMAND ${EXTRA_FASTQC_PARAMETERS} -t $PROCESSORS -o ./${OUTPUT}/ ${line}"

	scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${SAMPLE}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	chmodString="chmod 777 ${tempScript}"
	echo -e `${chmodString}`

	echo -e "source loadModules.sh\n\n" > ${tempScript}
	echo "$COMMAND_LINE" >> ${tempScript}

	if [ $SCHEDULER == "sge" ]; then
	    SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${fileName}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
	else
	    SUBMIT_COMMAND="bsub -q $QUEUE -J ${fileName}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}_${SAMPLE}_R${PAIREDNUMBER}.%J.out -e ${JOBS}/${COMMAND}_${SAMPLE}_R${PAIREDNUMBER}.%J.error bash ${tempScript}"
	fi

	date=`date`
	echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
	cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
	echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

	echo `${SUBMIT_COMMAND}`

    fi
done < "$findResults"  # Sends the $findResults file as input to the while-loop
