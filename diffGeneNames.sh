#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# diffExpXgeneName - custom script for Annotating CuffDiff gene_exp.diff files with Gene Names
#   - The script requires the file ensemblToGeneName.txt to run. It basically reads in the Ensemble id from the gene_exp.diff file,
#
#   - finds the gene name corresponding to the Ensemble id and appends the gene name to the respective line in the gene_exp.diff file.
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
PROCESSORS=1
MEMORY="4096"       # 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="72:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"        # short = max 8 hours;   long = max ? days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${GENENAMES_INPUT}"
OUTPUT=$(eval echo "$GENENAMES_OUTPUT")
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

#----GLOB INPUT----------------------------------------------------------------------------------------------------------
findResults=${OUTPUT}/input.list.txt
FIND_INPUT=$(eval echo "$GENENAMES_FIND_INPUT_DIR")
find ${FIND_INPUT} -name "${GENENAMES_FIND_INPUT_PATTERN}" > $findResults  # grabs all the "*-over-*" subdirs (cuffdiff output dirs)
#------------------------------------------------------------------------------------------------------------------------

COMMAND=perl

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a .diff pathName

    diffFileName=$(basename $line)  # strip all dirs
    diffFile=${line#*/}    # strip first dir
    overDir=$(basename $(dirname $line))

    scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${diffFileName}.XXXXXXXXXXX"
    echo -e "\n${scriptString}"
    tempScript=`${scriptString}`
    echo -e "\n${tempScript}"
    chmodString="chmod 777 ${tempScript}"
    echo -e `${chmodString}`


    COMMAND_LINE="${COMMAND} -f ${TAILOR}/diffExpXgeneName.pl ${TAILOR}/ensemblToGeneName.txt ${INPUT}/${diffFile} ${OUTPUT}/${diffFile}"

    echo -e "source loadModules.sh\n" > ${tempScript}
    echo -e "mkdir $OUTPUT/$overDir\n" >> ${tempScript}
    echo "$COMMAND_LINE" >> ${tempScript}

    if [ $SCHEDULER == "sge" ]; then
	SUBMIT_COMMAND="qsub ${EXTRA_SUBMIT_PARAMETERS} -q $QUEUE -cwd -S /bin/bash -N ${diffFileName}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
    else
	SUBMIT_COMMAND="bsub ${EXTRA_SUBMIT_PARAMETERS} -q $QUEUE -J ${diffFileName}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}.%J.out -e ${JOBS}/${COMMAND}.%J.error bash ${tempScript}"
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

