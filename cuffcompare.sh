#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# cuffcompare - help analyze the transfrags you assemble. The program cuffcompare helps you:
#
#    - Compare your assembled transcripts to a reference annotation
#    - Track Cufflinks transcripts across multiple samples (e.g. across a time course)
#
# This produces the following files in the "cuffcom-out" directory with the following files
#     cuffcmp.transcripts.gtf.tmap,
#     cuffcmp.transcripts.gtf.refmap,
#     cuffcmp.tracking,
#     cuffcmp.stats,
#     cuffcmp.loci,
#     cuffcmp.combined.gtf
#     cuffcmp.stats
#          Reports statistics related to the "accuracy" of the transcripts when compared to the reference annotation data.
#          Gene finding measures of "sensitivity" and "specificity" are calculated at various levels (nucleotide, exon, intron, transcript, gene)
#     cuffcmp.combined.gtf
#          Reports a GTF file containing the "union" of all transfrags in each sample
#
#  - Can send you emails when each job starts to execute and when each is finished so
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
MEMORY="4192"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="3:58"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="short"        # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT="${CUFFCOMPARE_INPUT}"
OUTPUT="${CUFFCOMPARE_OUTPUT}"
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
# find ${INPUT} -name "${CUFFCOMPARE_FIND_INPUT_PATTERN}" > ${OUTPUT}/assemblies.txt
#-----------------------------------------------------------------

COMMAND=cuffcompare

# COMMAND_LINE="${COMMAND} -r ${GTF} -R -V -s ${FASTADIR} -o $OUTPUT/cuffcmp -i ${OUTPUT}/assemblies.txt"
COMMAND_LINE="${COMMAND} ${EXTRA_CUFFCOMPARE_PARAMETERS} ../${CUFFCOMPARE_INPUT}/${CUFFCOMPARE_ASSEMBLIES_GTF}"

scriptString="mktemp -p ${SCRIPTS} ${COMMAND}.assemblies.XXXXXXXXXXX"
echo -e "\n${scriptString}"
tempScript=`${scriptString}`
echo -e "\n${tempScript}"
chmodString="chmod 777 ${tempScript}"
echo -e `${chmodString}`

echo -e "source loadModules.sh\n\n" > ${tempScript}
echo -e "cd ${OUTPUT}\n\n" >> ${tempScript}
echo "$COMMAND_LINE" >> ${tempScript}

if [ $SCHEDULER == "sge" ]; then
    SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${COMMAND} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
else
    SUBMIT_COMMAND="bsub -q $QUEUE -J ${COMMAND} -n ${PROCESSORS} -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -B -o ${JOBS}/${COMMAND}.%J.out -e ${JOBS}/${COMMAND}.%J.error bash ${tempScript}"
fi

date=`date`
echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

echo `${SUBMIT_COMMAND}`
# rm ${tempScript}
