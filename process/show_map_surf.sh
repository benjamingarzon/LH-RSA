#!/bin/sh
WD=/data/lv0/MotorSkill/fmriprep/analysis/higherlevel/Trained_Untrained/volume/
SUBJECTS_DIR=/data/lv0/MotorSkill/labels/subject
NAME=kk

for HEMI in rh lh; do
mri_vol2surf --mov $WD/$NAME.nii.gz \
   --sd $SUBJECTS_DIR \
   --regheader cvs_avg35_inMNI152 \
   --trgsubject cvs_avg35_inMNI152 \
   --hemi $HEMI \
   --surf white \
   --projfrac-avg 0.2 0.8 0.2 \
   --o $WD/$HEMI.$NAME.func.gii
done

freeview -f $SUBJECTS_DIR/cvs_avg35_inMNI152/surf/rh.inflated:overlay=$WD/rh.$NAME.func.gii \
$SUBJECTS_DIR/cvs_avg35_inMNI152/surf/lh.inflated:overlay=$WD/lh.$NAME.func.gii 
#   --surf-fwhm fwhm : smooth output surface (mm)

SUBJECT=lue1101
SESSION=3
run=3
HOMEDIR=/data/lv0/MotorSkill
FUNCREF=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-T1w_desc-preproc_bold.nii.gz #_space-T1w_boldref.nii.gz #
FUNCREFMEAN=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-T1w_desc-preproc_bold_mean.nii.gz
SUBJECTS_DIR=$HOMEDIR/fmriprep/freesurfer
fslmaths $FUNCREF -Tmean $FUNCREFMEAN
freeview $SUBJECTS_DIR/sub-$SUBJECT/mri/T1.mgz $FUNCREFMEAN -f $SUBJECTS_DIR/sub-$SUBJECT/surf/rh.pial $SUBJECTS_DIR/sub-$SUBJECT/surf/rh.white
mri_vol2surf --mov $WD/$NAME.nii.gz \
   --sd $SUBJECTS_DIR \
   --regheader cvs_avg35_inMNI152 \
   --trgsubject cvs_avg35_inMNI152 \
   --hemi $HEMI \
   --surf white \
   --projfrac-avg 0.2 0.8 0.2 \
   --o $WD/$HEMI.$NAME.func.gii
done

SUB=sub-lue5104
SESS=6
RUN=1
freeview fmriprep/$SUB/anat/${SUB}_desc-preproc_T1w.nii.gz fmriprep/$SUB/ses-${SESS}/func/${SUB}_ses-${SESS}_task-sequence_run-${RUN}_space-T1w_desc-preproc_bold.nii.gz 
#../data_BIDS/$SUB/ses-${SESS}/func/${SUB}_ses-${SESS}_task-sequence_run-0${RUN}_bold.nii.gz

