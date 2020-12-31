#!/bin/bash

#------ UMASS Boston - Riley Lab ------------------------------------------------------------------------------
# trim.sh - Trim adapters & nucleotides, performing Quality Control on the fastq files
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
MEMORY="1024"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="24:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"        # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-------------------
INPUT=$TRIM_INPUT
OUTPUT_FASTQC=${TRIM_OUTPUT_FASTQC}
OUTPUT_TRIM=${TRIM_OUTPUT_TRIM}
SCRIPTS=${OUTPUT_TRIM}/${JOBS_SCRIPTS_DIR}
JOBS=${OUTPUT_TRIM}/${JOBS_OUT_DIR}
#----------------------------

#----OUTPUT------------------
if [ ! -d ${OUTPUT_FASTQC} ]; then
    mkdir ${OUTPUT_FASTQC}
fi
if [ ! -d ${OUTPUT_TRIM} ]; then
    mkdir ${OUTPUT_TRIM}
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
cp "${TAILOR_CONFIG}" "${OUTPUT_TRIM}/${configFileName}.$(date +%F_%R)" 

#----GLOB INPUT--------------------------------
findResults=${OUTPUT_TRIM}/input.list.txt
find ${INPUT} -name "${TRIM_FIND_INPUT_PATTERN}" > $findResults
#----------------------------------------------

#----TRIM PARAMETERS--------
COMMAND=trim_galore
#---------------------------

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a fastq file

    fileName=$(basename $line)

    if [[ $fileName =~ $TRIM_INPUT_REGEX ]]; then
     if [[ "${BASH_REMATCH[2]}" == 1* ]]; then # starts with 1

	SAMPLE="${BASH_REMATCH[1]}"
	PAIREDNUMBER="${BASH_REMATCH[2]}"

	scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.${SAMPLE}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	chmodString="chmod 777 ${tempScript}"
	echo -e `${chmodString}`

        R1file=$line
        # R2file="${INPUT}/${SAMPLE}_R${PAIREDNUMBER}.fastq.gz"
        R2file="${INPUT}/${SAMPLE}_R2.fastq.gz"

# The following line was added for toxoplasmo analysis
#	R2file="${INPUT}/${SAMPLE}_R2.fastq"

	# echo "SAMPLE=${SAMPLE}; PAIREDNUMBER=${PAIREDNUMBER}; R1file=${R1file}; R2file=${R2file}"
	EXTRA_PARAMETERS=$(eval echo "$EXTRA_TRIM_PARAMETERS")
	echo -e "\n${EXTRA_PARAMETERS}"

    	COMMAND_LINE="$COMMAND $EXTRA_PARAMETERS --fastqc_args \"--outdir ${OUTPUT_FASTQC}/\" --output_dir ${OUTPUT_TRIM}"

	echo -e "source loadModules.sh\n\n" > ${tempScript}
	echo "$COMMAND_LINE" >> ${tempScript}

	if [ $SCHEDULER == "sge" ]; then
            SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${SAMPLE}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
	else
	    SUBMIT_COMMAND="bsub -q $QUEUE -J ${SAMPLE}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -B -o ${JOBS}/${COMMAND}_${SAMPLE}.%J.out -e ${JOBS}/${COMMAND}_${SAMPLE}.%J.error bash ${tempScript}"
	fi

	date=`date`
	echo -e "\n# $date\n"      >> ${OUTPUT_TRIM}/${COMMAND}.jobs.log
	echo -e "\n# Job Script\n" >> ${OUTPUT_TRIM}/${COMMAND}.jobs.log
	cat ${tempScript}          >> ${OUTPUT_TRIM}/${COMMAND}.jobs.log
	echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT_TRIM}/${COMMAND}.jobs.log
	echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT_TRIM}/${COMMAND}.jobs.log

	echo `${SUBMIT_COMMAND}`
	# rm ${tempScript}

    fi
   fi

done < "$findResults"  # Sends the $findResults file as input to the while-loop
