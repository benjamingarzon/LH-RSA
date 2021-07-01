#/bin/sh
# do resampling

SUBJECT=$1
SESS=$2
RUN1=1
RUN2=2
MYDIR=/home/xgarzb@GU.GU.SE/Data/LeftHand/Lund1/fmriprep/fmriprep/$SUBJECT/ses-$SESS/func/
cd $MYDIR

if []; then
mv ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz \
${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-brain_mask_orig.nii.gz 

mri_vol2vol --mov  ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-brain_mask_orig.nii.gz \
 --targ ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN2}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz \
 --o ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz  --regheader

mv ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-preproc_bold_orig.nii.gz 

mri_vol2vol --mov ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-preproc_bold_orig.nii.gz  \
 --targ ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN2}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
 --o ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz  --regheader

mv ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz \
${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold_orig.nii.gz

mri_vol2vol --mov ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold_orig.nii.gz \
 --targ ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN2}_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz \
 --o ${SUBJECT}_ses-${SESS}_task-sequence_run-${RUN1}_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz --regheader
fi

# for i in `find */sub-lue*/ses*/func/sub-lue*_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz`; do val=`fslval $i dim1`; if [ ! $val == 108 ]; then echo $i $val; fi; done
#./do_resampling_MNI.sh sub-lue2102 3 #
#./do_resampling_MNI.sh sub-lue2105 5 #
#./do_resampling.sh sub-lue2106 6
#./do_resampling.sh sub-lue3106 1
#./do_resampling.sh sub-lue4205 1
#./do_resampling.sh sub-lue5104 1
#./do_resampling_MNI.sh sub-lue3104 2
#./do_resampling_MNI.sh sub-lue3205 2
#./do_resampling_MNI.sh sub-lue3205 3


