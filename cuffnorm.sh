#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# cuffnorm - find significant changes in transcript expression, splicing, and promoter use
#
#   - Groups of replicates for each experiement (case) are defined in settings.sh
#
#   - Calls singleCuffnorm.sh for each pair of cases you want to compare
#
#   - Can Send you emails when each job starts to execute and when each is finished so
#     that you know when to submit the jobs for the next step
#
#   - NOTE: The easiest, least-error-prone way to create the groups below is to use the "ls -l" of the tophat_results
#           directory and then group the listed subdirectories accordingly
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

#----PATHS-----------------------------------------------------------
OUTPUT="${CUFFNORM_OUTPUT}"
#--------------------------------------------------------------------

#----OUTPUT-----------------
if [ ! -d ${OUTPUT} ]; then
    mkdir ${OUTPUT}
fi
#---------------------------

# Save this version of the config file
configFileName=$(basename ${TAILOR_CONFIG})
cp "${TAILOR_CONFIG}" "${OUTPUT}/${configFileName}.$(date +%F_%R)" 

#----------------------------------------------------------------------------------------------------------------------------------------
# Write the group arrays to file with entries separated by newlines
#----------------------------------------------------------------------------------------------------------------------------------------
for i in "${!experimentNames[@]}"; do

    EXPER_FILE="${OUTPUT}/group${i}.list.txt"

    experimentGroupString="experimentGroups${i}"
    eval experimentGroup=\( \${${experimentGroupString}[@]} \)

    printf '%s\n' "${experimentGroup[@]}"  > $EXPER_FILE

done


#----------------------------------------------------------------------------------------------------------------------------------------
# Call singleCuffNorm on each requested comparison
#
# NOTE: singleCuffnorm.sh sample1Group sample1Name sample2Group sample2Name - generates expression fold changes of sample2 over sample1.
#
# So for Case vs Control the call should be: "singleCuffnorm.sh controlGroup controlName caseGroup caseName".
#
#----------------------------------------------------------------------------------------------------------------------------------------
for (( c=0; c<$NUM_COMPARISONS; c++ )); do

    GROUP_1="cuffnorm${c}[0]"
    GROUP_2="cuffnorm${c}[1]"

    GROUP_1_VAR=${!GROUP_1}
    GROUP_2_VAR=${!GROUP_2}

    EXPER_1_NAME="${experimentNames[${GROUP_1_VAR}]}"
    EXPER_2_NAME="${experimentNames[${GROUP_2_VAR}]}"

    EXPER_1_FILE="${OUTPUT}/group${GROUP_1_VAR}.list.txt"
    EXPER_2_FILE="${OUTPUT}/group${GROUP_2_VAR}.list.txt"

    # EXPER_1_NAME="${experimentNames[{!GROUP_1}]}"
    # EXPER_2_NAME="${experimentNames[{!GROUP_2}]}"

    # EXPER_1_FILE="${OUTPUT}/group${!GROUP_1}.list.txt"
    # EXPER_2_FILE="${OUTPUT}/group${!GROUP_2}.list.txt"

    # echo -e "\n\nCuffnorm comparison $c"
    # echo -e "\n   GROUP_1_VAR=$GROUP_1_VAR, GROUP_2_VAR=$GROUP_2_VAR"
    # echo -e "\n   EXPER_1_NAME=$EXPER_1_NAME, EXPER_2_NAME=$EXPER_2_NAME"
    # echo -e "\n   EXPER_1_FILE=$EXPER_1_FILE, EXPER_2_FILE=$EXPER_2_FILE"

    singleCuffnorm.sh $EXPER_1_FILE $EXPER_1_NAME $EXPER_2_FILE $EXPER_2_NAME

done
