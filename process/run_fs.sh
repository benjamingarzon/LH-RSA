#!/bin/sh
# run freesurfers

EXPERT_FILE=/home/benjamin.garzon/Software/LeftHand/process/expert.opts

do_fs_cross(){
# check not done yet

SUB=$1
T1w=$2
PD=$3
T2w=$4

if [ ! -f "$SUBJECTS_DIR/$SUB/stats/lh.aparc.stats" ]; then 
    rm -r $SUBJECTS_DIR/$SUB/
    echo "$SUBJECTS_DIR/$SUB/ not done"
    
    recon-all -autorecon1 -s ${SUB}.PD -i $PD -hires -expert $EXPERT_FILE
    mri_vol2vol --mov $SUBJECTS_DIR/${SUB}.PD/mri/brainmask.mgz \
    --targ $SUBJECTS_DIR/${SUB}.PD/mri/rawavg.mgz --regheader --nearest \
    --o $SUBJECTS_DIR/${SUB}.PD/mri/brainmasknative.nii.gz

    fslmaths $SUBJECTS_DIR/${SUB}.PD/mri/brainmasknative.nii.gz -bin -ero $SUBJECTS_DIR/${SUB}.PD/mri/mask.nii.gz 
    fslmaths $T1w -mas $SUBJECTS_DIR/${SUB}.PD/mri/mask.nii.gz $SUBJECTS_DIR/${SUB}.PD/mri/T1wmasked.nii.gz 

    recon-all -autorecon1 -noskullstrip -s $SUB -i $SUBJECTS_DIR/${SUB}.PD/mri/T1wmasked.nii.gz -hires \
    -expert $EXPERT_FILE

    cd $SUBJECTS_DIR/$SUB/mri
    cp T1.mgz brainmask.auto.mgz
    ln -s brainmask.auto.mgz brainmask.mgz
    recon-all -autorecon2 -autorecon3 -s $SUB -T2 $T2w -T2pial -hires #-expert $EXPERT_FILE

    rm -r $SUBJECTS_DIR/${SUB}.PD 
    
else
    echo "$SUBJECTS_DIR/$SUB/ already done"
fi

}

export SUBJECTS_DIR=/home/share/MotorSkill/freesurfer
WD=/home/share/MotorSkill/data_BIDS/ #$1
SUBJECT=sub-lue001 #SUBJECT=$2 # add sub-

echo Doing subject $SUBJECT
ANATLIST=$WD/$SUBJECT/ses-*/anat

for DIR in $ANATLIST; do
    N=`echo $DIR | cut -d'/' -f8 | cut -d '-' -f2` 
    T1w=$DIR/MP2RAGE.nii.gz 
    PD=$DIR/mag.nii.gz 
    T2w=$DIR/${SUBJECT}_ses-${N}_T2w.nii.gz 
    
    if [ -e "$T1w" ] && [ -e "$PD" ] && [ -e "$T2w" ]; then
    
        # all available
        echo "Running cross-sectional"    
        do_fs_cross ${SUBJECT}.$N $T1w $PD $T2w
    
        # run longit
    #    do_fs_long ${SUBJECT} 
    
    fi

done