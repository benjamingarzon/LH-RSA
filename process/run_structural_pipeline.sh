#!/usr/bin/sh

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer
EXPERT_FILE=/home/benjamin.garzon/Software/LeftHand/process/expert.opts

# remove myelin inhomogeneity
# 1201.3, 1201.4, 1201.5 bad quality scan, not used for constructing base

MAXPROCS=15
NSESSIONS=7
SUBJECTS="lue3101 lue3102 lue3103 lue3104 lue3105 lue3106 lue3107 lue3201 lue3202 lue3203 lue3204 lue3205 lue3206 lue3207"
SUBJECTS="lue3205 lue3204 lue3203 lue3201"

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
VBM_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/vbm

OVERWRITE=0
DOFS=1
for SUBJECT in $SUBJECTS; do
    if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-*_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then    
       nice ./run_fs.sh $SUBJECTS_DIR $WD $VBM_DIR $SUBJECT $EXPERT_FILE $OVERWRITE $DOFS > $HOMEDIR/logs/${SUBJECT}.fs.log
    fi
    
    sleep 10
    NPROCS=`pgrep run_fs |wc -w`
    while [ $NPROCS -gt $MAXPROCS ]; do
      echo "$NPROCS processes running"
      sleep 1000
      NPROCS=`pgrep run_fs |wc -w`
    done
done

#WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
#rm $WD/*/ses*/anat/*template*  $WD/*/ses*/anat/*weights* $WD/*/ses*/anat/*.lta $WD/*/ses*/anat/*native* $WD/*/ses*/anat/*restore* $WD/*/ses*/anat/*masked*

