#!/usr/bin/sh

# missing data subjects: 
# lue1103 ses3 run 4 
# lue1104 ses-1 run 2-5



HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer

PHASES="2"

NSESSIONS=6 # wave 4
#NSESSIONS=7 

SESSIONS=`seq $NSESSIONS`
#SESSIONS="7 6 5 4 3 2"

SUBJECTS="lue4101 lue4102 lue4104 lue4106 lue4201 lue4202 lue4203 lue4204 lue4205 lue4206" # lue4207 lue4103 lue4105 lue4107"
#SUBJECTS="lue3101 lue3102 lue3103 lue3104 lue3105 lue3106 lue3107 lue3201 lue3202 lue3203 lue3204 lue3205 lue3206 lue3207"
#SUBJECTS="lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"


WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
RESPONSES_FILE=$HOMEDIR/responses/trials_table_Lund1fmri_clean.csv

NTRIALS=32
RESPONSES_OPTIONS="2.1 0.5 0 6.0 6"
CONFOUND_INDICES="26,27,28,29,30,31"
RUNS=0 # Set to 0 to  detect them from data

conda activate lhenv2

for SUBJECT in $SUBJECTS; do

if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-?_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then 

     for PHASE in $PHASES; do
       for SESSION in $SESSIONS; do
         echo  ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE "$RUNS" "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR 
         echo "Log file: $HOMEDIR/logs/${SUBJECT}.${SESSION}.phase${PHASE}.log"
         nice ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE $RUNS "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR > $HOMEDIR/logs/${SUBJECT}.${SESSION}.phase${PHASE}.log 
         if [ "$PHASE" -eq "3" ]; then
             # Do it only once
             break 1
         fi 
       done
     done

fi
done

