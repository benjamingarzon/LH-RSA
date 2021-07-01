#!/bin/sh
WD=/data/lv0/MotorSkill/fmriprep/analysis/higherlevel/Trained_Untrained/volume/tests/average
SUBJECTS_DIR=/data/lv0/MotorSkill/labels/subject
NAME=mean

if [ ]; then 
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


#$SUBJECTS_DIR=$HOMEDIR/fmriprep/fmriprep/
HEMI=rh
SUBJECT=lue1103
SESSION=3
run=3
HOMEDIR=/data/lv0/MotorSkill

FUNCREF=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz #_space-T1w_boldref.nii.gz #
FUNCREFMEAN=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-MNI152NLin2009cAsym_desc-preproc_bold_mean.nii.gz
#fslmaths $FUNCREF -Tmean $FUNCREFMEAN
for HEMI in rh lh; do
OUTFILE=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/$HEMI.test.func.gii
mri_vol2surf --mov $FUNCREF \
   --sd $SUBJECTS_DIR \
   --regheader cvs_avg35_inMNI152 \
   --trgsubject cvs_avg35_inMNI152 \
   --hemi $HEMI \
   --surf white \
   --projfrac-avg 0.2 0.8 0.2 \
   --o $OUTFILE

done

freeview $SUBJECTS_DIR/cvs_avg35_inMNI152/mri/T1.mgz $FUNCREF -f $SUBJECTS_DIR/cvs_avg35_inMNI152/surf/rh.pial:overlay=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/rh.test.func.gii $SUBJECTS_DIR/cvs_avg35_inMNI152/surf/lh.pial:overlay=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/lh.test.func.gii

exit 1
fi
HOMEDIR=/data/lv0/MotorSkill
SUBJECTS_DIR=$HOMEDIR/fmriprep/freesurfer
HEMI=rh
SUBJECT=lue1103
SESSION=3
run=3

FUNCREF=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-T1w_desc-preproc_bold.nii.gz #_space-T1w_boldref.nii.gz #
FUNCREFMEAN=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-T1w_desc-preproc_bold.nii.gz
#fslmaths $FUNCREF -Tmean $FUNCREFMEAN


if [ 1 ]; then
for HEMI in rh lh; do
mri_vol2surf --mov $FUNCREF \
   --sd $SUBJECTS_DIR \
   --regheader sub-$SUBJECT \
   --trgsubject sub-$SUBJECT  \
   --hemi $HEMI \
   --surf white \
   --projfrac-avg 0.2 0.7 0.1 \
   --o $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-individual-hemi-${HEMI}_bold.func.gii

mri_surf2surf --sd $SUBJECTS_DIR \
   --srcsubject sub-${SUBJECT} \
   --srcsurfval $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-individual-hemi-${HEMI}_bold.func.gii \
   --trgsubject fsaverage6\
   --srcsurfreg sphere.reg \
   --trgsurfval $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-fsaverage6-hemi-${HEMI}_bold.func.gii \
   --hemi $HEMI \
   --cortex
done
fi

freeview $SUBJECTS_DIR/sub-$SUBJECT/mri/T1.mgz $FUNCREF -f $SUBJECTS_DIR/sub-$SUBJECT/surf/rh.pial:overlay=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-individual-hemi-rh_bold.func.gii $SUBJECTS_DIR/sub-$SUBJECT/surf/rh.white

freeview -f /usr/local/freesurfer/7.1.1-1/subjects/fsaverage6/surf/rh.inflated:overlay=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-fsaverage6_hemi-R_bold.func.gii:overlay=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-${run}_space-fsaverage6-hemi-${HEMI}_bold.func.gii 
exit 1 
SUB=sub-lue5104
SESS=6
RUN=1
freeview fmriprep/$SUB/anat/${SUB}_desc-preproc_T1w.nii.gz fmriprep/$SUB/ses-${SESS}/func/${SUB}_ses-${SESS}_task-sequence_run-${RUN}_space-T1w_desc-preproc_bold.nii.gz 
#../data_BIDS/$SUB/ses-${SESS}/func/${SUB}_ses-${SESS}_task-sequence_run-0${RUN}_bold.nii.gz

