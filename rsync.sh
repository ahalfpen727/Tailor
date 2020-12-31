#!/bin/bash

#---------------------------------------------------------------------------------------
# RSYNC all files from SOURCE machine to the local OUTPUT directory
# WARNING - YOU MUST CHANGE $SOURCE_USER TO YOU AND CREATE AND CREATE THE PASSWORD FILE!
#  - Note: requires the sshpass utility
#  - can send you emails when each job starts to execute and when each is finished so that you know when to submit the jobs for the next step
# This script can also be used to download Genome builds, Fasta files, reads
# from any repository in a more effective manner than ftp
#---------------------------------------------------------------------------------------

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
PROCESSORS=1
MEMORY="4096"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="192:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"        # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----LOCAL PATHS-------------
OUTPUT=$RSYNC_OUTPUT
SCRIPTS=${OUTPUT}/${JOBS_SCRIPTS_DIR}
JOBS=${OUTPUT}/${JOBS_OUT_DIR}
#---------------------------

#----Output-----------------
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

COMMAND=rsync

#----Glob Input on the remote host using ssh---------------------------------------------------------------
findResults=${OUTPUT}/input.list.txt
sshpass -f ${SOURCE_PASS_FILE} ssh ${SOURCE_USER}@{SOURCE_HOST} "find ${SOURCE_DIR} -name '$SOURCE_TYPE'" > $findResults
#----------------------------------------------------------------------------------------------------------

while read filename; do

    COMMAND_LINE="${COMMAND} ${EXTRA_RSYNC_PARAMETERS} --password-file=${SOURCE_PASS_FILE} ${SOURCE_USER}@{SOURCE_HOST}:${SOURCE_DIR}/${filename} $OUTPUT"

    scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${filename}.XXXXXXXXXXX"
    echo -e "\n${scriptString}"
    tempScript=`${scriptString}`
    echo -e "\n${tempScript}"
    chmodString="chmod 777 ${tempScript}"
    echo -e `${chmodString}`

    echo -e "source loadModules.sh\n\n" > ${tempScript}
    echo "$COMMAND_LINE" >> ${tempScript}

    if [ $SCHEDULER == "sge" ]; then
	SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
    else
	SUBMIT_COMMAND="bsub -q $QUEUE -J ${command} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${command}.%J.out -e ${JOBS}/${command}.%J.error bash ${tempScript}"
    fi

    date=`date`
    echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
    cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

    echo `${SUBMIT_COMMAND}`
    # rm ${tempScript}

done < "$findResults"  # Sends the $findResults file as input to the while-loop
