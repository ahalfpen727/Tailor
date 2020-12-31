#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
#  launchR.eachDiffFile.sh - launches R scripts for each cuffdiff *.diff file in each cuffdiff output directory
#
#   - $1 must contain the prefix of the R script
#
#   - Can send you emails when each job starts to execute and when each is finished so
#     that you know when to submit the jobs for the next step
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
MEMORY="60192"       # 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="3:58"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"        # short = max 8 hours;   long = max ? days                                                                                                                                             
#QUEUE="short"        # short = max 8 hours;   long = max ? days
#export QUEUE
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
# INPUT="${GENENAMES}${INPUT_LABEL}"
INPUT="${LAUNCHR_INPUT}"
OUTPUT=$(eval echo "$LAUNCHR_OUTPUT")
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
FIND_INPUT=$(eval echo "$LAUNCHR_EACHDIFFFILE_FIND_INPUT_DIR")
find ${FIND_INPUT} -name "${LAUNCHR_EACHDIFFFILE_FIND_INPUT_PATTERN}" > $findResults
# grabs all the "*-over-*" subdirs (cuffdiff output dirs)
#------------------------------------------------------------------------------------------------------------------------

COMMAND=R

echo -e "\nINPUT=$INPUT; OUTPUT=$OUTPUT; findResults=$findResults";

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a cuffdiff output directory

    diffFile=$(basename $line)
    diffFileLength=${#diffFile}
    diffDir=${diffFile:0:$(($diffFileLength - 5))} # strip off the ".diff" extension
    overDir=$(basename $(dirname $line))
    # diffDir=$(basename $line)

    [[ $overDir =~ $LAUNCHR_INPUT_REGEX ]]
    over=${BASH_REMATCH[1]}
    under=${BASH_REMATCH[2]}

    echo -e "\ndiffFile=${diffFile}; diffDir=${diffDir}; overDir=${overDir}; over=${over}; under=${under}"

    scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${diffDir}.XXXXXXXXXXX"
    echo -e "\n${scriptString}"
    tempScript=`${scriptString}`
    echo -e "\n${tempScript}"
    chmodString="chmod 777 ${tempScript}"
    echo -e `${chmodString}`

    if [ ! -d "${OUTPUT}/${overDir}/${diffDir}" ]; then
	mkdir "${OUTPUT}"
	mkdir "${OUTPUT}/${overDir}"
	mkdir "${OUTPUT}/${overDir}/${diffDir}"
    fi

    # COMMAND_LINE="${XVFB_RUN} ${COMMAND} CMD BATCH --no-site-file --no-init-file --no-restore --slave ${1}.R ${INPUT}/${overDir}/${diffFile} ${OUTPUT}/${overDir}/${diffDir}"
    # COMMAND_LINE="${XVFB_RUN} ${COMMAND} CMD BATCH --no-site-file --no-init-file --no-restore --slave '--args inFile=\"${diffFile}\" inDir=\"${INPUT}/${overDir}\" outDir=\"${OUTPUT}/${diffDir}\" over=\"${over}\" under=\"${under}\"' ${1}.R ${OUTPUT}/${diffDir}/${1}.Rout"

#export GO_DIFF_EXPRESSED_ALPHA_VALUE
#export GO_HYPER_GEO_ALPHA_VALUE

    COMMAND_LINE="${XVFB_RUN} ${COMMAND} CMD BATCH --no-restore --no-save '--args inFile=\"${diffFile}\" inDir=\"${INPUT}/${overDir}\" outDir=\"${OUTPUT}/${overDir}/${diffDir}\" over=\"${over}\" under=\"${under}\"' ${TAILOR}/${1}.R ${OUTPUT}/${overDir}/${diffDir}/${1}.Rout"

    echo -e "env\n\n" > ${tempScript}
    echo -e "source loadModules.sh\n\n" >> ${tempScript}
    echo "$COMMAND_LINE" >> ${tempScript}

    if [ $SCHEDULER == "sge" ]; then
	SUBMIT_COMMAND="qsub ${EXTRA_SUBMIT_PARAMETERS} -q $QUEUE -cwd -S /bin/bash -N ${diffDir}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
    else
	SUBMIT_COMMAND="bsub ${EXTRA_SUBMIT_PARAMETERS} -q $QUEUE -J ${diffDir}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}.%J.out -e ${JOBS}/${COMMAND}.%J.error bash ${tempScript}"
    fi

    date=`date`
    echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
    cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

    echo `${SUBMIT_COMMAND}`
    # rm ${tempScript}

    # exit

done < "$findResults"  # Sends the $findResults file as input to the while-loop

