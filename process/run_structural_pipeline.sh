#!/usr/bin/sh

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer
EXPERT_FILE=/home/benjamin.garzon/Software/LeftHand/process/expert.opts

# remove myelin inhomogeneit
# T2 do manual correg
# correct 1206 neck, maybe multiply by T2??

NSESSIONS=7
SUBJECTS="lue1101 lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
SUBJECTS="lue1101 lue1103 lue1104"
WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
VBM_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/vbm

OVERWRITE=1
DOFS=1
for SUBJECT in $SUBJECTS; do
    if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-*_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then    
       echo ./run_fs.sh $SUBJECTS_DIR $WD $SUBJECT $EXPERT_FILE
       nice ./run_fs.sh $SUBJECTS_DIR $WD $VBM_DIR $SUBJECT $EXPERT_FILE $OVERWRITE > $HOMEDIR/logs/${SUBJECT}.fs.log &
    fi
done

#WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
#rm $WD/*/ses*/anat/*template*  $WD/*/ses*/anat/*weights* $WD/*/ses*/anat/*.lta $WD/*/ses*/anat/*native* $WD/*/ses*/anat/*restore* $WD/*/ses*/anat/*masked*