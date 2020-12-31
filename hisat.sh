#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# HISAT2 - insert aware aligner serving as the replacement to tophat
# aligns reads to the genome and calculates alignments, junctions, insertions, and deletions
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
#PROCESSORS=6
PROCESSORS=9
MEMORY="7592"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
#MEMORY="5592"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="3:59"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
#DURATION="96:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"        # short = max 4 hours;   long = max 30 days
#QUEUE="long"        # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${HISAT_INPUT}"
OUTPUT="${HISAT_OUTPUT}"
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
if [ ! -d ${ANNOTATION_INDEX} ]; then
    mkdir ${ANNOTATION_INDEX}
fi
if [ ! -d ${OUTPUT}/${SAMPLE_OUT} ]; then
    mkdir ${OUTPUT}/${SAMPLE_OUT}
fi
#---------------------------

# Save this version of the config file
configFileName=$(basename ${TAILOR_CONFIG})
cp "${TAILOR_CONFIG}" "${OUTPUT}/${configFileName}.$(date +%F_%R)" 

#----GLOB INPUT--------------------------------
findResults=${OUTPUT}/input.list.txt
find ${INPUT} -name "${HISAT_FIND_INPUT_PATTERN}" > $findResults
#----------------------------------------------

COMMAND=hisat2

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a fastq file

    fileName=$(basename $line)

    if [[ $fileName =~ $HISAT_INPUT_REGEX ]]; then

	i=${BASH_REMATCH[1]}

	# Cannonical Command
	COMMAND_LINE="cd ${OUTPUT};  touch ${OUTPUT}/${i}_out_mapped.sam;  ${COMMAND} -p ${PROCESSORS} ${EXTRA_PARAMETERS} -1 ${INPUT}/${i}_R1_val_1.fq.gz -2 ${INPUT}/${i}_R2_val_2.fq.gz  -S ${OUTPUT}/${i}_out_mapped.sam"


	scriptString="mktemp -p ${SCRIPTS}  ${COMMAND}.${i}.XXXXXXXXXXX"
	echo -e "\n${scriptString}"
	tempScript=`${scriptString}`
	echo -e "\n${tempScript}"
	chmodString="chmod 777 ${tempScript}"
	echo -e `${chmodString}`

	EXTRA_PARAMETERS=$(eval echo "$EXTRA_HISAT_PARAMETERS")
	echo -e "\n${EXTRA_PARAMETERS}"


   ### When running Tophat, supply a reference genome annotation file (GTF or GFF).  The transcriptome index and bowtie2 index can then be built ##
   ### with the -G option pointing to the GTF or GFF genome annotation and  the --transcriptome-index option pointing to the directory that will ##
   ### house the newly built transcriptome.  This initial tophat run can be perfomed without performing mapping/aligning to reduce CPU time. ######

	echo -e "source loadModules.sh\n\n" > ${tempScript}
	echo "$COMMAND_LINE" >> ${tempScript}

#	echo -e "grep -i multiple ./${OUTPUT}/${SAMPLE}_out/align_summary.txt >> ${OUTPUT}/multiple_mapping.txt\n" >> ${tempScript}
#	echo -e "grep -i concordant ./${OUTPUT}/${SAMPLE}_out/align_summary.txt >> ${OUTPUT}/concordant_mapping.txt\n" >> ${tempScript}
#	echo -e "grep -i overall ./${OUTPUT}/${SAMPLE}_out/align_summary.txt >> ${OUTPUT}/overall_mapping.txt\n" >> ${tempScript}

	if [ $SCHEDULER == "sge" ]; then
            SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${i}-${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
	else
	    SUBMIT_COMMAND="bsub -q $QUEUE -J ${i}-${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -B -o ${JOBS}/${COMMAND}-${i}.%J.out -e ${JOBS}/${COMMAND}-${i}.%J.error bash ${tempScript}"
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
