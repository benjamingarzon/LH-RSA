#!/usr/bin/sh


PHASE=2
HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer
#SUBJECTS="lue1101 lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
SUBJECTS="lue1201 lue1204"
SUBJECTS="lue1105"

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
RESPONSES_FILE=$HOMEDIR/responses/trials_table_Lund1fmri_clean.csv
NSESSIONS=7
NTRIALS=32
RESPONSES_OPTIONS="2.1 0.5 0 6.0 6"
CONFOUND_INDICES="26,27,28,29,30,31"
RUNS=0 # Set to 0 to  detect them from data



if [ ]; then

for SUBJECT in $SUBJECTS; do
 if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-?_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then 

   if [ "$PHASE" -eq "3" ]; then
     ./pipeline.sh $SUBJECT 1 $RESPONSES_FILE $NTRIALS $PHASE $RUNS "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR > $HOMEDIR/logs/${SUBJECT}.phase3.log
   else
     for SESSION in `seq $NSESSIONS`; do
       echo  ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE $RUNS "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR
       ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE $RUNS "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR > $HOMEDIR/logs/${SUBJECT}.${SESSION}.phase${PHASE}.log 
     done
   fi
 fi
done

fi

#exit 1

for SUBJECT in $SUBJECTS; do
    if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-*_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then    
       echo ./run_fs.sh $SUBJECTS_DIR $WD $SUBJECT
       ./run_fs.sh $SUBJECTS_DIR $WD $SUBJECT > $HOMEDIR/logs/${SUBJECT}.fs.log &
    fi
done
exit 1

