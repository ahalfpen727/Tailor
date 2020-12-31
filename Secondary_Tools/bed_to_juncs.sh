#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# bed_to_juncs 
#  
#  a tool in the BEDtools packages can be used to analyze tophat output
#  This script helps create 5 bp window around splice site junctions for input into IGV
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
MEMORY="8912"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="94:00"  # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"       # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${BEDJUNCS_INPUT}"
OUTPUT="${BEDJUNCS_OUTPUT}"
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
 find ${INPUT} -type d -maxdepth 1 -name "${BEDJUNCS_FIND_INPUT_PATTERN}" > $findResults  # search for only directories and no subdirs
#----------------------------------------------

COMMAND=bedtools

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a directory that contains a accepted_hit.bam file
    fileName=$(basename $line)

    if [[ $fileName =~ $BEDJUNCS_INPUT_REGEX ]]; then

	i=${BASH_REMATCH[1]}

#	COMMAND_LINE="${COMMAND} -p ${PROCESSORS} ${EXTRA_CUFFLINKS_PARAMETERS} -o ${OUTPUT}/${i}_out ${INPUT}/${i}_out/accepted_hits.bam"
	COMMAND_LINE="cat ${INPUT}/${i}_out/junctions.bed | ${COMMAND} bed_to_juncs | ${COMMAND} flankBed -g ${GENOME_SIZE_FILE} -b 5 > ${OUTPUT}/${i}_out/junction_windows.txt"

	scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${i}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	chmodString="chmod 777 ${tempScript}"
	echo -e `${chmodString}`

	echo -e "source loadModules.sh\n\n" > ${tempScript}
	echo "$COMMAND_LINE" >> ${tempScript}

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
