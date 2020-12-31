#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# cuffmerge - allows you to merge together several Cufflinks assemblies.
# It also handles running Cuffcompare for you, and automatically filters a number of transfrags that are
# probably artfifacts. If you have a reference GTF file available, you can provide it to the script in order to
# gracefully merge novel isoforms and known isoforms and maximize overall assembly quality. The main purpose of
# this script is to make it easier to make an assembly GTF file suitable for use with Cuffdiff.
#
#   - AUTOMATICALLY CALLS CUFFCOMPARE
#
#   - Sends you emails when each job starts to execute and when each is finished so
#    that you know when to submit the jobs for the next step
#
#   - YOU MUST SET THE $LSB_MAILTO TO YOUR EMAIL ADDRESS!!
#        example: export LSB_MAILTO=todd.riley@umb.edu
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
MEMORY="4096"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="3:58"  # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"       # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${CUFFMERGE_INPUT}"
OUTPUT="${CUFFMERGE_OUTPUT}"
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

#----GLOB INPUT----------------------------------------------------
 find ${INPUT} -name "${CUFFMERGE_FIND_INPUT_PATTERN}" > ${OUTPUT}/${CUFFMERGE_FIND_OUTPUT_FILE}
#------------------------------------------------------------------

COMMAND=cuffmerge

COMMAND_LINE="${COMMAND} ${EXTRA_CUFFMERGE_PARAMETERS} -p ${PROCESSORS} -o $OUTPUT ${OUTPUT}/${CUFFMERGE_FIND_OUTPUT_FILE}"

scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.assemblies.XXXXXXXXXXX"
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
    SUBMIT_COMMAND="bsub -q $QUEUE -J ${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}.%J.out -e ${JOBS}/${COMMAND}.%J.error bash ${tempScript}"
fi

date=`date`
echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

echo `${SUBMIT_COMMAND}`
# rm ${tempScript}
