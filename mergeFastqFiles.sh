#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# MergeFastqFiles.sh - merge fastqfiles into 1 for each sample
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
PROCESSORS=1
MEMORY="8192"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="192:00"  # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"       # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT=${MERGE_INPUT}
OUTPUT=${MERGE_OUTPUT}
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

cp $findResults "${findResults}.old"
rm -f $findResults
for DIR in $(find ${INPUT} -type d -name "Project_*" -not \( -path "*/Temp/*" -prune \)); do \
    echo -e "${DIR}\n"; \
    find $DIR -type d -name "${MERGE_FIND_INPUT_PATTERN}" -print  >> $findResults; \
done
# find ${INPUT} -type d -maxdepth 1 -name "*" > $findResults  # search for only directories and no subdirs
#----------------------------------------------

COMMAND=zcat

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a directory that contains a accepted_hit.bam file

    files=$(ls $line)
    filesArray=($files)
    line1=${filesArray[0]}
    # echo -e "\n\nline = $line\n\n"
    echo -e "\n\nline1 = $line1\n\n"
    # exit 1

    if [[ $line1 =~ ${MERGE_REGEX} ]]; then

	dir=$(basename ${line})
	# dir=${BASH_REMATCH[1]}
	sample=${BASH_REMATCH[$SAMPLE_INDEX]}
	barcode=${BASH_REMATCH[$BARCODE_INDEX]}
	lane=${BASH_REMATCH[${LANE_INDEX}]}
	# read=${BASH_REMATCH[4]}
	# split=${BASH_REMATCH[5]}

	echo -e "\nline1 = $line1"
	echo -e "\ndir = $dir"
	echo -e "\nsample = $sample"
	echo -e "\nbarcode = $barcode"
	echo -e "\nlane = $lane"
	# echo -e "\nread = $read"
	# echo -e "\nsplit = $split"
        # exit 1

	# eval function evaluates the string which will perform variable replacement
	EVAL_R1_GLOB=$(eval echo $R1_GLOB)
	EVAL_R2_GLOB=$(eval echo $R2_GLOB)
	EVAL_R1_DESTINATION=$(eval echo $R1_DESTINATION)
	EVAL_R2_DESTINATION=$(eval echo $R2_DESTINATION)

	COMMAND_LINE1="${COMMAND} ${EXTRA_MERGE_PARAMETERS} ${line}/${EVAL_R1_GLOB}  | gzip -c > ${OUTPUT}/${EVAL_R1_DESTINATION}"
	COMMAND_LINE2="${COMMAND} ${EXTRA_MERGE_PARAMETERS} ${line}/${EVAL_R2_GLOB}  | gzip -c > ${OUTPUT}/${EVAL_R2_DESTINATION}"

	scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${i}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	chmodString="chmod 777 ${tempScript}"
	echo -e `${chmodString}`

	echo -e "source loadModules.sh\n\n" > ${tempScript}
	echo "$COMMAND_LINE1" >> ${tempScript}
	echo "$COMMAND_LINE2" >> ${tempScript}

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
