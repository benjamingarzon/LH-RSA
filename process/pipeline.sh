# fix fieldmaps
# align runs well
# use mrqc instead of mriprep to get good segmentation
# fmriprep single subject?
# try without using derivatives in the regressors...
# how long should the trials be?

#STEPS
# Create EVs from response file
# Convert MRI data to BIDS
# Process with fmriprep to native space
# Create fsl setup files (remove derivatives?)
# Process fMRI in native space
# Open copes and organize them 
# RSA


HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lundpilot1
WD=$HOMEDIR/fmriprep
WORK=$HOMEDIR/work
ECHOTIME1=0.001544
ECHOTIME2=0.002544 
# + 0.0038523

s=102
export NTRIALS=27
ORIG_FSF_FILE=$HOMEDIR/fMRInoreg_multistretch.fsf

s=101
export NTRIALS=24
ORIG_FSF_FILE=$HOMEDIR/fMRInoreg_multi.fsf

###############################################
# Create explanatory variables for fMRI analysis
###############################################
cd $HOMEDIR/responses/sub-$s
rm -r run*
#python ~/Software/LeftHand/process/organize_responses.py trialsfile-lup1s001_fmri.csv 2.8 0.5 0 23 4
python ~/Software/LeftHand/process/organize_responses.py trialsfile-lup0s002_fmri_edited.csv 2.8 0.5 0 23 4


###############################################
# Convert to BIDS
###############################################
source activate lhenv

cd $HOMEDIR/

heudiconv -d '7T033/7T033{subject}_*/*/Dicom/*/*.dcm' -s $s -f /home/benjamin.garzon/Software/LeftHand/process/heuristic_conv.py -o data_BIDS_vol/sub-$s -b --overwrite
mv data_BIDS_vol/sub-$s/fmap/sub-${s}_epi1.nii.gz data_BIDS_vol/sub-$s/fmap/sub-${s}_magnitude.nii.gz 

fslmaths data_BIDS_vol/sub-$s/fmap/sub-${s}_epi2.nii.gz -mul 6.28 data_BIDS_vol/sub-$s/fmap/sub-${s}_fieldmap_rads.nii.gz 
fugue --loadfmap=data_BIDS_vol/sub-$s/fmap/sub-${s}_fieldmap_rads.nii.gz -m --savefmap=data_BIDS_vol/sub-$s/fmap/sub-${s}_fieldmap.nii.gz

for run in 1 2 3 4; do
sed -i 's/Axis/Direction/g' func/sub-${s}_task-sequence_run-0${run}_bold.json
done 

echo -e "{\n\"PhaseEncodingDirection\": \"j\", \n\"Units\": \"rad/s\", \n\"IntendedFor\": [\"func/sub-${s}_task-sequence_run-01_bold.nii.gz\", \"func/sub-${s}_task-sequence_run-02_bold.nii.gz\", \"func/sub-${s}_task-sequence_run-03_bold.nii.gz\", \"func/sub-${s}_task-sequence_run-04_bold.nii.gz\"]\n}" > data_BIDS_vol/sub-$s/fmap/sub-${s}_fieldmap.json 

#echo -e "{\n\"EchoTime1\": $ECHOTIME1, \n\"EchoTime2\": $ECHOTIME2, \n\"IntendedFor\": [\"func/sub-${s}_task-sequence_run-01_bold.nii.gz\", \"func/sub-${s}_task-sequence_run-02_bold.nii.gz\", \"func/sub-${s}_task-sequence_run-03_bold.nii.gz\", \"func/sub-${s}_task-sequence_run-04_bold.nii.gz\"]\n}" > data_BIDS_vol/sub-$s/fmap/sub-${s}_phasediff.json 

rm -r data_BIDS_vol/sourcedata


###############################################
# Create singularity file container
###############################################

#cd /tmp
#docker run --privileged -t --rm \
#    -v /var/run/docker.sock:/var/run/docker.sock \
#    -v /tmp/fmriprep:/output \
#    singularityware/docker2singularity \
#    poldracklab/fmriprep:latest


#cp /tmp/fmriprep/poldracklab_fmriprep_latest-*.img $HOMEDIR/fmriprep.img

#rm -r $WD
rm -r $WORK
mkdir $WD
mkdir $WORK
cp /usr/local/freesurfer/license.txt $WD/



###############################################
# Process with fmriprep
###############################################
PYTHONPATH="" singularity run \
   $HOMEDIR/fmriprep.img \
   $HOMEDIR/data_BIDS_vol \
   $WD \
   participant \
   --ignore slicetiming \
   --ignore fieldmaps \
   --fs-license-file \
   $WD/license.txt \
   --nthreads 4 \
   -w $WORK \
   --write-graph \
   --output-space 'T1w' \
   --participant-label $s 

#   --use-syn-sdc \
#   --use-aroma \
#   --output-space 'fsnative' \

#for s in $SUBJECTS;  
#do
#  PYTHONPATH="" singularity run $HOMEDIR/processed/fmriprep.img $HOMEDIR/data_BIDS $WD participant --use-aroma --use-syn-sdc --fs-license-file $WD/license.txt --nthreads 4 -w $WORK --participant-label $s & 
#done

# collect the data
fslmerge -t $HOMEDIR/fmriprep/analysis/sub-${s}/data.nii.gz $HOMEDIR/fmriprep/fmriprep/sub-${s}/func/sub-${s}_task-sequence_run-*_bold_space-T1w_preproc.nii.gz
fslmaths $HOMEDIR/fmriprep/analysis/sub-${s}/data.nii.gz -Tmean $HOMEDIR/fmriprep/analysis/sub-${s}/data_mean.nii.gz


###############################################
# Run fMRI analysis
###############################################

# using raw (no fmriprep)
STRUCT=$HOMEDIR/data_BIDS_vol/sub-${s}/anat/sub-${s}_T1w_brain.nii.gz

NEVS=3
for run in 2 3 4 1; 
do
  FUNC=$HOMEDIR/data_BIDS_vol/sub-${s}/func/sub-${s}_task-sequence_run-0${run}_bold.nii.gz
  RUN_DIR=$HOMEDIR/data_BIDS_vol/sub-${s}/func/run$run/
  rm -r $RUN_DIR
  mkdir $RUN_DIR

  FSF_FILE=$RUN_DIR/fMRI.fsf
  
  cp $HOMEDIR/fMRI.fsf $FSF_FILE
  cp $HOMEDIR/responses/run$run/*.csv $RUN_DIR
  
  cd $RUN_DIR

  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@STRUCTURAL_BRAIN%$STRUCT%" $FSF_FILE
  sed -i "s%@ANALYSIS%$RUN_DIR/analysis%" $FSF_FILE  
  
  for EV in `seq $NEVS`; 
  do  
     sed -i "s%@EV$EV%$RUN_DIR/EV$EV.csv%" $FSF_FILE        
  done

  # now run it
  feat $FSF_FILE&

done


export NTRIALS=24
# remove dummies!
cd $HOMEDIR/responses
rm -r run*
python ~/Software/LeftHand/process/organize_responses.py trialsfile-lup0s002_fmri_edited.csv 2.8 0.5 23 4
cd $HOMEDIR/
~/Software/LeftHand/process/setup_GLM.sh


# single trial analysis
STRUCT=$HOMEDIR/data_BIDS_vol/sub-${s}/anat/sub-${s}_T1w_brain.nii.gz

for run in 1 2 3 4; 
do
  FUNC=$HOMEDIR/data_BIDS_vol/sub-${s}/func/sub-${s}_task-sequence_run-0${run}_bold.nii.gz
  RUN_DIR=$HOMEDIR/data_BIDS_vol/sub-${s}/func/run_single$run/
  rm -r $RUN_DIR
  mkdir $RUN_DIR

  FSF_FILE=$RUN_DIR/fMRI.fsf
  
  cp $HOMEDIR/fMRIsingle_complete.fsf $FSF_FILE
  cp $HOMEDIR/responses/run$run/SINGLE* $RUN_DIR
  
  cd $RUN_DIR

  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@STRUCTURAL_BRAIN%$STRUCT%" $FSF_FILE
  sed -i "s%@ANALYSIS%$RUN_DIR/analysis%" $FSF_FILE  
  sed -i "s%@NEVS%$(($NTRIALS + 1))%" $FSF_FILE
  sed -i "s%@DOUBLENEVS%$((2*$NTRIALS + 2))%" $FSF_FILE
  
  for EV in `seq $(($NTRIALS + 1))`; 
  do  
     sed -i "s%@EV$EV%$RUN_DIR/SINGLE${EV}%" $FSF_FILE        
  done

  # now run it
  feat $FSF_FILE&

done


# gather PE files in 1

# using fmriprep outputs
SUBDIR=$HOMEDIR/fmriprep/analysis/sub-${s}

cd $SUBDIR
for run in 1 2 3 4; 
do

  FUNC=$HOMEDIR/fmriprep/fmriprep/sub-${s}/func/sub-${s}_task-sequence_run-0${run}_bold_space-T1w_preproc.nii.gz

  RUN_DIR=$SUBDIR/run$run
  rm -r $RUN_DIR
  mkdir $RUN_DIR 

    # get confounds 
  cut -d$'\t' -f 28,29,30,31,32,33 $HOMEDIR/fmriprep/fmriprep/sub-${s}/func/sub-${s}_task-sequence_run-0${run}_bold_confounds.tsv | tail -n +2 > $RUN_DIR/mc_run-0${run}.csv
  
  for TRIAL in `seq $NTRIALS`; 
#  for TRIAL in 1 2; 
#  for TRIAL in `seq 24 27`; 
 
   do
   echo $TRIAL

  MODEL_DIR=$RUN_DIR/model$TRIAL
  mkdir $MODEL_DIR

  FSF_FILE=$MODEL_DIR/fMRI.fsf
  
  cp $ORIG_FSF_FILE $FSF_FILE
  cp $HOMEDIR/responses/sub-${s}/run$run/TRIAL$TRIAL $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${s}/run$run/OTHER$TRIAL $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${s}/run$run/FIXATION $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${s}/run$run/STRETCH $MODEL_DIR/
  
  cd $MODEL_DIR

  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@ANALYSIS%$MODEL_DIR/analysis%" $FSF_FILE
  
  sed -i "s%@EV1%$MODEL_DIR/TRIAL$TRIAL%" $FSF_FILE        
  sed -i "s%@EV2%$MODEL_DIR/OTHER$TRIAL%" $FSF_FILE        
  sed -i "s%@EV3%$MODEL_DIR/FIXATION%" $FSF_FILE        
  sed -i "s%@EV4%$MODEL_DIR/STRETCH%" $FSF_FILE        
  sed -i "s%@CONFOUNDS%$RUN_DIR/mc_run-0${run}.csv%" $FSF_FILE        

  # now run it
  feat $FSF_FILE 
  done &
done

# optimize High-pass? xx

FX_FILE=tstat
# merge parameters
cd $HOMEDIR/fmriprep/analysis/sub-${s}
for run in 1 2 3 4; do
  fslmerge -t run$run/effects_$run `ls -v run$run/model*/analysis.feat/stats/${FX_FILE}1.nii.gz` 
  fslmerge -t run$run/derivatives_$run `ls -v run$run/model*/analysis.feat/stats/${FX_FILE}2.nii.gz` 
done
  fslmerge -t effects `ls -v run*/effects*` 
  fslmerge -t derivatives `ls -v run*/derivatives*` 

  fslmerge -t other `ls -v run$run/model*/analysis.feat/stats/pe3.nii.gz` 
  fslmerge -t fixation `ls -v run$run/model*/analysis.feat/stats/pe5.nii.gz` 
  fslmerge -t stretch `ls -v run$run/model*/analysis.feat/stats/pe7.nii.gz` 
  
  fslmaths other -Tmean -s 3 other_mean
  fslmaths fixation -Tmean -s 3 fixation_mean
  fslmaths effects -Tmean -s 3 effects_mean
  fslmaths derivatives -Tmean -s 3 derivatives_mean
  fslmaths stretch -Tmean -s 3 stretch_mean
  
# prewhiten the parameters 
#fslmaths res4d.nii.gz -sqr -Tmean -mul `fslnvols res4d.nii.gz` -div `cat dof` -sqrt res_sigma
fslmaths sigmasquared.nii.gz -sqrt sigma.nii.gz


###############################################
# Run searchlight analysis
###############################################

# prepare surfaces
#export PATH=$PATH:/usr/local/AFNI
#pymvpa2-prep-afni-surf -d $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf -e $HOMEDIR/fmriprep/analysis/sub-${s}/effects.nii.gz -r  $HOMEDIR/fmriprep/analysis/sub-${s}/suma

RADIUS=10.0
NPROC=35
mri_vol2vol --mov $HOMEDIR/fmriprep/fmriprep/sub-${s}/anat/sub-${s}_T1w_brainmask.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${s}/effects.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${s}/mask.nii.gz --nearest
fslmaths $HOMEDIR/fmriprep/analysis/sub-${s}/effects.nii.gz -abs -Tmin -thr 0 -mas $HOMEDIR/fmriprep/analysis/sub-${s}/mask -bin $HOMEDIR/fmriprep/analysis/sub-${s}/mask

mkdir $HOMEDIR/fmriprep/analysis/sub-${s}/surf
# to facilitate visualization
mris_convert --to-scanner $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/lh.pial $HOMEDIR/fmriprep/analysis/sub-${s}/surf/lh.pial.gii
mris_convert --to-scanner $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/rh.pial $HOMEDIR/fmriprep/analysis/sub-${s}/surf/rh.pial.gii
mris_convert --to-scanner $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/lh.white $HOMEDIR/fmriprep/analysis/sub-${s}/surf/lh.white.gii
mris_convert --to-scanner $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/rh.white $HOMEDIR/fmriprep/analysis/sub-${s}/surf/rh.white.gii
mris_convert --to-scanner $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/lh.midthickness $HOMEDIR/fmriprep/analysis/sub-${s}/surf/lh.midthickness.gii
mris_convert --to-scanner $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/rh.midthickness $HOMEDIR/fmriprep/analysis/sub-${s}/surf/rh.midthickness.gii

ln -s $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/lh.inflated $HOMEDIR/fmriprep/analysis/sub-${s}/surf/lh.inflated
ln -s $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/rh.inflated $HOMEDIR/fmriprep/analysis/sub-${s}/surf/rh.inflated
ln -s $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/lh.curv $HOMEDIR/fmriprep/analysis/sub-${s}/surf/lh.curv
ln -s $HOMEDIR/fmriprep/freesurfer/sub-${s}/surf/rh.curv $HOMEDIR/fmriprep/analysis/sub-${s}/surf/rh.curv
ln -s $HOMEDIR/fmriprep/freesurfer/fsaverage/surf/lh.inflated $HOMEDIR/fmriprep/analysis/sub-${s}/surf/lh.fsaverage.inflated
ln -s $HOMEDIR/fmriprep/freesurfer/fsaverage/surf/rh.inflated $HOMEDIR/fmriprep/analysis/sub-${s}/surf/rh.fsaverage.inflated
ln -s $HOMEDIR/fmriprep/freesurfer/fsaverage/surf/lh.curv $HOMEDIR/fmriprep/analysis/sub-${s}/surf/lh.fsaverage.curv
ln -s $HOMEDIR/fmriprep/freesurfer/fsaverage/surf/rh.curv $HOMEDIR/fmriprep/analysis/sub-${s}/surf/rh.fsaverage.curv


SUBJECTS_DIR=$HOMEDIR/fmriprep/freesurfer

# prepare engine
for hemi in rh lh; do
python ~/Software/LeftHand/process/prepare_queryengine.py \
    $HOMEDIR/fmriprep/analysis/sub-${s} \
    $HOMEDIR/responses/sub-${s}/sequences.csv \
    $hemi $RADIUS $NPROC
done


# run searchlight
for hemi in rh lh; do
python ~/Software/LeftHand/process/surface_searchlight.py \
    $HOMEDIR/fmriprep/analysis/sub-${s}  \
    $HOMEDIR/responses/sub-${s}/sequences.csv \
    $hemi $RADIUS $NPROC

#for hemi in lh rh; do
# resample to average
mri_convert $SUBDIR/surf/$hemi.sl_accuracy_${RADIUS}.func.gii $SUBDIR/surf/$hemi.sl_accuracy_${RADIUS}.func.mgh
mri_surf2surf --srcsubject sub-${s} \
    --srcsurfval $SUBDIR/surf/$hemi.sl_accuracy_${RADIUS}.func.mgh \
    --trgsubject fsaverage \
    --trgsurfval $SUBDIR/surf/$hemi.sl_accuracy_${RADIUS}.fsaverage.func.mgh \
    --hemi $hemi \
    --fwhm-trg 5 \
    --cortex

mri_convert $SUBDIR/surf/$hemi.sl_accuracy_${RADIUS}.fsaverage.func.mgh $SUBDIR/surf/$hemi.sl_accuracy_${RADIUS}.fsaverage.func.gii

done
rm $SUBDIR/surf/*.mgh

# visualize
freeview -f rh.inflated:overlay=rh.sl_accuracy_10.0.func.gii lh.inflated:overlay=lh.sl_accuracy_10.0.func.gii
freeview -f rh.fsaverage.inflated:overlay=rh.sl_accuracy_10.0.fsaverage.func.gii:curv=rh.fsaverage.curv \
lh.fsaverage.inflated:overlay=lh.sl_accuracy_10.0.fsaverage.func.gii:curv=lh.fsaverage.curv

mkdir $SUBDIR/results

###############################################
# Extract roi data
###############################################



