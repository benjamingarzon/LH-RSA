#!/usr/bin/sh

HOMEDIR=/data/lv0/MotorSkill #~/Data/LeftHand/Lund1
SUBJECTS_DIR=/data/lv0/MotorSkill/freesurfer #~/Data/LeftHand/Lund1/freesurfer

export NPROC=20
export MAXFEATPROCS=20
export FMRIPREPPROCS=20


PHASES="8"

NSESSIONS=7

SESSIONS=`seq $NSESSIONS`

SUBJECTS="lue5101 lue5102 lue5103 lue5104 lue5105 lue5106 lue5107 lue5201 lue5202 lue5203 lue5204 lue5205 lue5206 lue5207"
#SUBJECTS="lue3101 lue3102 lue3103 lue3104 lue3105 lue3106 lue3107 lue3201 lue3202 lue3203 lue3204 lue3205 lue3206 lue3207"
SUBJECTS="lue4101 lue4102 lue4103 lue4104 lue4105 lue4106 lue4107 lue4201 lue4202 lue4203 lue4204 lue4205 lue4206 lue4207"
#SUBJECTS="lue2101 lue2102 lue2103 lue2104 lue2105 lue2106 lue2107 lue2201 lue2202 lue2203 lue2204 lue2205 lue2206 lue2207"

#SUBJECTS="lue1101 lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
#SUBJECTS="lue1201"
SUBJECTS="lue1101"

CWD=`pwd`
WD=~/Data/LeftHand/Lund1/data_BIDS

if [ "$SUBJECTS" ]; then
echo "Doing $SUBJECTS"
else

cd $WD 
SUBJECTS=`echo sub-lue* | sed 's/sub-//g'`
echo "Doing $SUBJECTS"
cd $CWD
fi

RESPONSES_FILE=$HOMEDIR/responses/complete_trials_fMRI_table.csv #trials_table_Lund1fmri_clean.csv

NTRIALS=32
RESPONSES_OPTIONS="2.1 0.5 0 6.0 6"
CONFOUND_INDICES="26,27,28,29,30,31"
RUNS=0 # Set to 0 to  detect them from data

conda init bash
conda activate lhenv2

for SUBJECT in $SUBJECTS; do

#if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-?_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then 

     for PHASE in $PHASES; do
       for SESSION in $SESSIONS; do
         if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/anat/sub-${SUBJECT}_ses-?_MP2RAGE.nii.gz | wc -l` -eq 1 ]; then 
           echo  ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE "$RUNS" "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR 
           echo "Log file: $HOMEDIR/logs/${SUBJECT}.${SESSION}.phase${PHASE}.log"
            nice ./pipeline.sh $SUBJECT $SESSION $RESPONSES_FILE $NTRIALS $PHASE $RUNS "$RESPONSES_OPTIONS" "$CONFOUND_INDICES" $HOMEDIR > $HOMEDIR/logs/${SUBJECT}.${SESSION}.phase${PHASE}.log 
            if [ "$PHASE" -eq "3" ]; then
               # Do it only once
               break 1
            fi
         fi 
       done
     done
#fi
done

