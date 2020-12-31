#!/bin/bash

#---------------------------------------------------------------------------------------------------
# Create hardlinks to the original fastq files. The purpose of the hardlinks is to replace long,
# meaningless file names wth short, meaningful ones - which aids in troubleshooting and validation.
#
#  - The links should have VERY SHORT meaningful names constructed from the following keywords:
#
#       CTRL = control
#       CASE = case
#       EXP  = experiment
#       SAMP = sample
#       SYMP = symptomatic
#       LN   = lane
#       MPLX = multiplex
#
#  - Can send you emails when each job starts to execute and when each is finished so
#    that you know when to submit the jobs for the next step
#
#---------------------------------------------------------------------------------------------------

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

#----PATHS------------------------------------------------------------------------------------------
INPUT=$LINKS_INPUT
OUTPUT=$LINKS_OUTPUT
#---------------------------------------------------------------------------------------------------

#----OUTPUT-----------------
if [ ! -d ${OUTPUT} ]; then
    mkdir ${OUTPUT}
fi
#---------------------------

# Save this version of the config file
configFileName=$(basename ${TAILOR_CONFIG})
cp "${TAILOR_CONFIG}" "${OUTPUT}/${configFileName}.$(date +%F_%R)" 

idx=0
while read -a linkMatrix$idx; do    # each line is a "source<tab>destination" fastq file pairing
    let idx++;
done < "$LINKS_FILE"  # tab separated table where 1st column is the source file and the 2nd column is the destination link


for (( c=0; c<$idx; c++ )); do

    sourceFile="linkMatrix${c}[0]"
    destinationFile="linkMatrix${c}[1]"

    ln ${EXTRA_LINK_PARAMETERS} ${INPUT}/${!sourceFile} ${OUTPUT}/${!destinationFile}

done


