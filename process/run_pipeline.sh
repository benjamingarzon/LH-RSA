#!/usr/bin/sh

# missing data subjects: 
# lue1103 ses3 run 4 
# lue1104 ses-1 run 2-5

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer

PHASES="2"
PHASE=2

NSESSIONS=7
SESSIONS=`seq $NSESSIONS`
SESSIONS="1"

#SUBJECTS="lue1107"
#SUBJECTS="lue1101"
#SUBJECTS="lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
SUBJECTS="lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
SUBJECTS="lue1202"

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
RESPONSES_FILE=$HOMEDIR/responses/trials_table_Lund1fmri_clean.csv

NTRIALS=32
RESPONSES_OPTIONS="2.1 0.5 0 6.0 6"
CONFOUND_INDICES="26,27,28,29,30,31"
RUNS=0 # Set to 0 to  detect them from data

conda activate lhenv2


for SUBJECT in $SUBJECTS; do

if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-?_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then 

   if [ "$PHASE" -eq "3" ]; then
     nice ./pipeline.sh $SUBJECT 1 $RESPONSES_FILE $NTRIALS $PHASE $RUNS "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR > $HOMEDIR/logs/${SUBJECT}.phase3.log
   else
     for SESSION in $SESSIONS; do
       for PHASE in $PHASES; do
         echo  ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE "$RUNS" "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR 
         echo "Log file: $HOMEDIR/logs/${SUBJECT}.${SESSION}.phase${PHASE}.log"
         nice ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE $RUNS "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR > $HOMEDIR/logs/${SUBJECT}.${SESSION}.phase${PHASE}.log 
       done
     done
   fi
 fi
done

