#!/bin/bash
#--------------------------------------------------------------------------------------------------------------
# BAMTOOLS stats on the accepted hits bam files not available in SAMTOOLS
 # (like coverage)
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
PROCESSORS=5
MEMORY="4096"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="72:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"        # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${BAMTOOLS_INPUT}"
OUTPUT="${BAMTOOLS_OUTPUT}"
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

#----GLOB INPUT---------------------------------------------------
findResults=${OUTPUT}/input.list.txt
find ${INPUT} -type d -maxdepth 1 -name "${BAMTOOLS_FIND_INPUT_PATTERN}" > $findResults  # search for only directories and no subdirs

# find ${INPUT} -name "accepted_hits.bam"  > $findResults

#-----------------------------------------------------------------

COMMAND=bamtools

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a directory that contains a accepted_hit.bam file

    if [[ $line =~ $BAMTOOLS_INPUT_REGEX ]]; then

        SAMPLE=$(basename ${BASH_REMATCH[1]})
        # SAMPLE=($basename ${BASH_REMATCH[1]})

        scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${SAMPLE}.XXXXXXXXXXX"
        echo -e "\n${scriptString}"
        tempScript=`${scriptString}`
        echo -e "\n${tempScript}"
        chmodString="chmod 777 ${tempScript}"
        echo -e `${chmodString}`

        COMMAND_LINE="${COMMAND} ${BAMTOOLS_COMMAND} ${EXTRA_BAMTOOLS_PARAMETERS} ${INPUT}/${SAMPLE}_out/accepted_hits.bam > ${OUTPUT}/${SAMPLE}_out.${BAMTOOLS_COMMAND}"


# ${BAMTOOLS_COMMAND} ${EXTRA_BAMTOOLS_PARAMETERS} $BAMTOOLS_INPUT_REGEX 
# ${BAMTOOLS_FIND_INPUT_PATTERN} "${BAMTOOLS_INPUT}" "${BAMTOOLS_OUTPUT}"



        echo -e "source ./loadModules.sh\n\n" > ${tempScript}
        echo "$COMMAND_LINE" >> ${tempScript}

        if [ $SCHEDULER == "sge" ]; then
             SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${SAMPLE}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas  ${tempScript}"
         else
             SUBMIT_COMMAND="bsub -q $QUEUE -J ${SAMPLE}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -B -o ${JOBS}/${COMMAND}.%J.out -e ${JOBS}/${COMMAND}.%J.error bash ${tempScript}"
        fi

        date=`date`
        echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
        echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
        cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
        echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
        echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

        echo `${SUBMIT_COMMAND}`
        # rm ${tempScript}


    fi


done < "$findResults"  # Sends the $findResults file as input to the while-loop
