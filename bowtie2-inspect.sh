#!/bin/bash

#--BOWTIE2-INSPECT-----------------------------------------
# pass in 'build' or 'inspect' as 1st word of the variable BOWTIE_EXTRA_PARAMETERS
# <bt2_base> The basename of the index files to write. By default, bowtie2-build writes files named NAME.1.bt2, NAME.2.bt2, NAME.3.# bt2, NAME.4.bt2, NAME.rev.1.bt2, and NAME.rev.2.bt2, where NAME is <bt2_base>.
# bowtie2-inspect
# This command takes in a bowtie2-made index and outputs the reference sequences used to create it.
#-----------------------------------------------------------------------------
#----SETTINGS---------------
source "${TAILOR_CONFIG}"
#---------------------------

#----JOB SUBMISSION PARAMETERS---------------------------------------------------------------------
PROCESSORS=4
MEMORY="3096"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="3:55"  # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"       # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${BOWTIE_INPUT}"
OUTPUT="${BOWTIE_OUTPUT}"
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

COMMAND=bowtie2-inspect

COMMAND_LINE="${COMMAND} ${INPUT} > ${OUTPUT}"

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

