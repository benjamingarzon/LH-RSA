#!/usr/bin/sh

# missing data subjects: 
# lue1103 ses3 run 4 
# lue1104 ses-1 run 2-5

# effects not done
# lue1106
#lue1201 3-7
#lue1202 
#lue1203 5-7

# redo 1207.1 run2 with fmriprep and run
# redo 1202 fmriprep
# after redone 1203 fmriprep and run

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer

#PHASES="8 9 10"
PHASES="4 5 6 7 8 9 10 11 12"
PHASES="1 4 5 6 7 8 9 10 11 12"
PHASES="1 4 5 6 7"
#PHASES="1 4 5"
PHASES="3"
#PHASES="5"
#PHASES="9 10 11"
#PHASES="10 11"


NSESSIONS=7
SESSIONS=`seq $NSESSIONS`
#SESSIONS="7 6 5 4 3 2 1"

# check lue1106 1201 1202 1203
SUBJECTS="lue1101 lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
#SUBJECTS="lue1101 lue1103 lue1104 lue1105 lue1106 lue1201"
#SUBJECTS="lue1207 lue1206 lue1205 lue1204 lue1203 lue1202"
#SUBJECTS="lue1106 lue1201 lue1204" 
SUBJECTS="lue1207 lue1202" # lue1204" #lue1202 lue1203 
#SUBJECTS="lue1203" # lue1204" #lue1202 lue1203 
#SUBJECTS="lue1106"

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

