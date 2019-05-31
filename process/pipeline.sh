#/bin/bash

# fix fieldmaps
# align runs well
# use mrqc instead of mriprep to get good segmentation
# try without using derivatives in the regressors...
# how long should the trials be?

#STEPS
# Create EVs from response file
# Convert MRI data to BIDS
# Process with fmriprep to native space
# Create fsl setup files
# Process fMRI in native space
# Open copes and organize them 
# Searchlight analyses
#CHANGE ANALYSIS OUTPUT DIR??


# ADAPT number of frames
# remove 'incorrect' from fMRI analysis
# detect runs properly; IMPORTANT TO SYNC THE RUNS WITH THE BEHAVIOUR!
# use expert option in Freesurfer?
# check MP2RAGE mask, 10 % needed?
# clean up == remove run data
# simplify what we are doing 

# PHASE = 0 do everything
###############################################
# Variable definitions
###############################################
SUBJECT=$1
SESSION=$2
RESPONSES_FILE=$3 #trialsfile-lup2s008_fmri.csv
export NTRIALS=$4 #28 # $2
PHASE=$5
RUNS=$6
ORGANIZE_FILE=$7
RESPONSES_OPTIONS=$8
CONFOUND_INDICES=$9 #28,29,30,31,32,33

HOMEDIR=/home/share/MotorSkill
#RESPONSES_OPTIONS="2.8 0.5 0 25 6" old version
#RESPONSES_OPTIONS="2.1 0.5 0 6.0"
PROGDIR=~/Software/LeftHand/process/
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

# Same always
SIMPLE3_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_3stretch.fsf
SIMPLE4_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_4stretch.fsf
MULTI_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_multistretch.fsf

WD=$HOMEDIR/fmriprep
SUBJECTS_DIR=$WD/freesurfer
WORK=$HOMEDIR/work_sub-${SUBJECT}_ses-${SESSION}
ECHOTIME1=0.001544
ECHOTIME2=0.002544 
TE=25 #ms
EPIFAC=39
WFS=22.366
SENSEP=3
FIELDSTRENGTH=7
EES=`echo "((1000 * $WFS)/($FIELDSTRENGTH*3.4*42.57 * ($EPIFAC+1))/$SENSEP)" | bc -l | awk '{printf "%f", $0}'`

INTENDEDFOR="[\"func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-01_bold.nii.gz\", \"func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-02_bold.nii.gz\", \"func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-03_bold.nii.gz\", \"func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-04_bold.nii.gz\"]"

FX_FILE=tstat
RADIUS=15.0
NPROC=35
NEVS=4
FWHM=5
SIGMA=`echo $FWHM/2.3548 | bc -l`

#s=101
#export NTRIALS=24
#ORIG_FSF_FILE=$HOMEDIR/fMRInoreg_multi.fsf

conda activate lhenv2
if [ ! -e $WD ]; then 
    mkdir $WD
fi


###############################################
# Create explanatory variables for fMRI analysis
###############################################
if [ $PHASE == 0 ] || [ $PHASE == 1 ]; then

cd $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}
rm -r run*
python $PROGDIR/$ORGANIZE_FILE $RESPONSES_FILE $RESPONSES_OPTIONS > CorrectTrials.csv

fi


###############################################
# Prepare fieldmaps
###############################################
if [ $PHASE == 0 ] || [ $PHASE == 2 ]; then

cd $HOMEDIR/
cp data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_epi1.nii.gz data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_magnitude.nii.gz 

fslmaths data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_epi2.nii.gz -mul 6.28 data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap_rads.nii.gz 
fugue --loadfmap=data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap_rads.nii.gz -m --savefmap=data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap.nii.gz

for run in $RUNS; do
    sed -i 's/Axis/Direction/g' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json
    myline=`cat data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json | grep EchoTrainLength`
    sed -i '/'"$myline"'/a \  \"EffectiveEchoSpacing\": '"$EES"',' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json    
done 

echo -e "{\n\"EchoTime1\": $ECHOTIME1, \n\"EchoTime2\": $ECHOTIME2, \n\"PhaseEncodingDirection\": \"j\", \n\"Units\": \"rad/s\", \n\"IntendedFor\": $INTENDEDFOR\n}" > data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap.json 

cd data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/
#fslroi sub-${SUBJECT}_ses-${SESSION}_T1w1.nii.gz re.nii.gz 0 1 
#fslroi sub-${SUBJECT}_ses-${SESSION}_T1w1.nii.gz im.nii.gz 1 1 
ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec2.nii.gz re.nii.gz 
ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.nii.gz im.nii.gz

fslroi sub-${SUBJECT}_ses-${SESSION}_MP2RAGE.nii.gz mag 1 1
#fslroi sub-${SUBJECT}_ses-${SESSION}_MP2RAGE.nii.gz phase 0 1

cd $PROGDIR/mprageconvert/
matlab -nosplash -nodisplay -r "addpath $PROGDIR/mprageconvert/Nifti_tools; create_mp2rage_command('$HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/', 're.nii.gz ', 'im.nii.gz', 'sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.json'); exit"

cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/
#fslmaths MP2RAGE.nii.gz -add 0.5 -mul 1000 -mas mask -uthr 999 sub-${SUBJECT}_ses-${SESSION}_T1w.nii.gz -odt int
cp MP2RAGE.nii.gz sub-${SUBJECT}_ses-${SESSION}_T1w.nii.gz
#fslmaths MP2RAGE.nii.gz -add 0.5 -mul mag sub-${SUBJECT}_ses-${SESSION}_PDT1w.nii.gz -odt int

#mv phase.nii.gz sub-${SUBJECT}_ses-${SESSION}_T2w.nii.gz
#cp sub-${SUBJECT}_ses-${SESSION}_T1w.json sub-${SUBJECT}_ses-${SESSION}_T2w.json
#fslmaths MP2RAGE.nii.gz -add 0.5 -mul 1000 -mul mag mp2mod
#recon-all -s recon -i mp2mod.nii.gz -autorecon1
#recon-all -s recon_m -i MP2RAGE.nii.gz -autorecon1
# add T2
fi
###############################################
# Create singularity file container
###############################################

if [ $PHASE == 0 ] || [ $PHASE == 3 ]; then

if [ ! -e $HOMEDIR/fmriprep.img ]; then
    cd /tmp
    docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /tmp/fmriprep:/output \
        singularityware/docker2singularity \
        poldracklab/fmriprep:latest
    
    cp /tmp/fmriprep/poldracklab_fmriprep_latest-*.img $HOMEDIR/fmriprep.img
fi



###############################################
# Process with fmriprep
###############################################


rm -r $WORK
mkdir $WORK
cp /usr/local/freesurfer/license.txt $WD/

#SUBJECTS_DIR=
#recon-all -s recon -i sub-106_ses-1_T1w.nii.gz -T2 sub-106_ses-1_T1w.nii.gz -T2pial -all

PYTHONPATH="" singularity run \
   $HOMEDIR/fmriprep.img \
   $HOMEDIR/data_BIDS \
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
   --participant-label ${SUBJECT} 

PYTHONPATH="" singularity run \
   $HOMEDIR/fmriprep.img \
   $HOMEDIR/data_BIDS \
   $WD \
   participant \
   --ignore slicetiming \
   --ignore fieldmaps \
   --fs-license-file \
   $WD/license.txt \
   --nthreads 4 \
   -w $WORK \
   --write-graph \
   --output-space 'fsnative' \
   --participant-label ${SUBJECT}   

rm -r $WORK

fi # PHASE 3
   
# --ignore fieldmaps \
#   --use-syn-sdc \
#   --use-aroma \


###############################################
# Run fMRI analysis
###############################################
# Native surface space from fmriprep
#STRUCT=$HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_T1w_brain.nii.gz

SUBDIR=$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/

if [ $PHASE == 0 ] || [ $PHASE == 4 ]; then

if [ ! -e $SUBDIR ]; then
mkdir $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/
mkdir $SUBDIR
mkdir $SUBDIR/surf

# collect the data
fslmerge -t $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data.nii.gz $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-*_bold_space-T1w_preproc.nii.gz
fslmaths $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data.nii.gz -Tmean $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data_mean.nii.gz

## correct session
# To facilitate visualization and projection to cortical surface
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.pial $SUBDIR/surf/lh.pial.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.pial $SUBDIR/surf/rh.pial.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.white $SUBDIR/surf/lh.white.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.white $SUBDIR/surf/rh.white.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.midthickness $SUBDIR/surf/lh.midthickness.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.midthickness $SUBDIR/surf/rh.midthickness.gii
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.sphere.reg $SUBDIR/surf/lh.sphere.reg
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.sphere.reg $SUBDIR/surf/rh.sphere.reg

fi

fi #PHASE 4


# Compute activation maps
if [ $PHASE == 0 ] || [ $PHASE == 5 ]; then

cd $SUBDIR

for xxx in 1; do # just to make it wait
for run in $RUNS; 
do

  FUNCVOL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_space-T1w_preproc.nii.gz
  FUNCSURFL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_space-fsnative.L.func.gii
  FUNCSURFR=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_space-fsnative.R.func.gii

  RUN_DIR=$SUBDIR/run$run

  rm -r $RUN_DIR
  mkdir $RUN_DIR

  #Smoothing 
  fslmaths $FUNCVOL -s $FWHM $RUN_DIR/smoothed_data.nii.gz
  #mri_surf2surf --s sub-${SUBJECT} --sval $FUNCSURFL --hemi lh --fwhm-trg $FWHM --tval $RUN_DIR/lh.smoothed.func.gii
  #mri_surf2surf --s sub-${SUBJECT} --sval $FUNCSURFR --hemi rh --fwhm-trg $FWHM --tval $RUN_DIR/rh.smoothed.func.gii
#  $HCPDIR/wb_command -metric-dilate $RUN_DIR/lh.smoothed.func.gii $SUBDIR/surf/lh.midthickness.gii 100 $RUN_DIR/lh.smoothed.func.gii -nearest 
#  $HCPDIR/wb_command -metric-dilate $RUN_DIR/rh.smoothed.func.gii $SUBDIR/surf/rh.midthickness.gii 100 $RUN_DIR/rh.smoothed.func.gii -nearest 

  cut -d$'\t' -f $CONFOUND_INDICES $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_confounds.tsv | tail -n +2 > $RUN_DIR/mc_run-0${run}.csv

  FSF_FILE=$RUN_DIR/fMRI.fsf
  
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/*.csv $MODEL_DIR/$RUN_DIR/
  
  cd $RUN_DIR

  if [ `cat EV2.csv | wc -l` == 0 ]; then
   cp $SIMPLE3_FSF_FILE $FSF_FILE
   sed -i "s%@EV1%$RUN_DIR/EV1.csv%" $FSF_FILE        
   sed -i "s%@EV3%$RUN_DIR/EV3.csv%" $FSF_FILE        
   sed -i "s%@EV4%$RUN_DIR/EV4.csv%" $FSF_FILE        
  else
   cp $SIMPLE4_FSF_FILE $FSF_FILE
   sed -i "s%@EV1%$RUN_DIR/EV1.csv%" $FSF_FILE        
   sed -i "s%@EV2%$RUN_DIR/EV2.csv%" $FSF_FILE        
   sed -i "s%@EV3%$RUN_DIR/EV3.csv%" $FSF_FILE        
   sed -i "s%@EV4%$RUN_DIR/EV4.csv%" $FSF_FILE          
  fi
  
  NVOLS=`fslnvols $RUN_DIR/smoothed_data.nii.gz`
  sed -i "s%@NVOLS%$NVOLS%" $FSF_FILE
  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@STRUCTURAL_BRAIN%$STRUCT%" $FSF_FILE
  sed -i "s%@ANALYSIS%$RUN_DIR/analysis%" $FSF_FILE  
  sed -i "s%@CONFOUNDS%$RUN_DIR/mc_run-0${run}.csv%" $FSF_FILE        
  
  # now run it
  feat_model fMRI
  rm -r $RUN_DIR/surfL $RUN_DIR/surfR $RUN_DIR/volume 
#  film_gls --rn=$RUN_DIR/surfL --sa --in=$RUN_DIR/lh.smoothed.func.gii --pd=fMRI.mat --con=fMRI.con --mode=surface --in2=$SUBDIR/surf/lh.midthickness.gii
#  film_gls --rn=$RUN_DIR/surfR --sa --in=$RUN_DIR/rh.smoothed.func.gii --pd=fMRI.mat --con=fMRI.con --mode=surface --in2=$SUBDIR/surf/rh.midthickness.gii 
  film_gls --rn=$RUN_DIR/volume --sa --ms=$FWHM --in=$RUN_DIR/smoothed_data.nii.gz --thr=1000 --pd=fMRI.mat --con=fMRI.con --mode=volume &
  
done

done
fi # PHASE 5

  ##### do second level analysis


# Prepare parameter estimates. Using fmriprep outputs
# assume PHASE5 already run
if [ $PHASE == 0 ] || [ $PHASE == 6 ]; then

for run in $RUNS; 
do

  FUNC=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_space-T1w_preproc.nii.gz

  RUN_DIR=$SUBDIR/run$run
#  rm -r $RUN_DIR
#  mkdir $RUN_DIR 

    # Get confounds 
  for TRIAL in `seq $NTRIALS`; 
  do
  echo $TRIAL

  MODEL_DIR=$RUN_DIR/model$TRIAL
  rm -r $MODEL_DIR
  mkdir $MODEL_DIR
  FSF_FILE=$MODEL_DIR/fMRI.fsf
  
  cp $MULTI_FSF_FILE $FSF_FILE
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/TRIAL$TRIAL $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/OTHER?_$TRIAL $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/FIXATION $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/STRETCH $MODEL_DIR/
  
  cd $MODEL_DIR
  NVOLS=`fslnvols $FUNC`
  sed -i "s%@NVOLS%$NVOLS%" $FSF_FILE
  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@ANALYSIS%$MODEL_DIR/analysis%" $FSF_FILE 
  sed -i "s%@EV1%$MODEL_DIR/TRIAL$TRIAL%" $FSF_FILE        
  sed -i "s%@EV2%$MODEL_DIR/OTHER1_$TRIAL%" $FSF_FILE        
  sed -i "s%@EV3%$MODEL_DIR/OTHER2_$TRIAL%" $FSF_FILE        
  sed -i "s%@EV4%$MODEL_DIR/OTHER3_$TRIAL%" $FSF_FILE        
  sed -i "s%@EV5%$MODEL_DIR/OTHER4_$TRIAL%" $FSF_FILE        
  sed -i "s%@EV6%$MODEL_DIR/FIXATION%" $FSF_FILE        
  sed -i "s%@EV7%$MODEL_DIR/STRETCH%" $FSF_FILE        
  sed -i "s%@CONFOUNDS%$RUN_DIR/mc_run-0${run}.csv%" $FSF_FILE        

  # Now run it
  feat $FSF_FILE 
  done &
done

fi # PHASE 6

if [ $PHASE == 0 ] || [ $PHASE == 7 ]; then

# Merge parameters
cd $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/

for run in $RUNS; do
  fslmerge -t run$run/effects_$run `ls -v run$run/model*/analysis.feat/stats/${FX_FILE}1.nii.gz` 
  fslmerge -t run$run/derivatives_$run `ls -v run$run/model*/analysis.feat/stats/${FX_FILE}2.nii.gz` 
done
  fslmerge -t effects `ls -v run*/effects*` 
  fslmerge -t derivatives `ls -v run*/derivatives*` 

#  fslmerge -t other1 `ls -v run$run/model*/analysis.feat/stats/pe3.nii.gz` 
#  fslmerge -t other2 `ls -v run$run/model*/analysis.feat/stats/pe5.nii.gz` 
#  fslmerge -t other3 `ls -v run$run/model*/analysis.feat/stats/pe7.nii.gz` 
#  fslmerge -t other4 `ls -v run$run/model*/analysis.feat/stats/pe9.nii.gz` 

  fslmerge -t fixation `ls -v run$run/model*/analysis.feat/stats/pe11.nii.gz` 
  fslmerge -t stretch `ls -v run$run/model*/analysis.feat/stats/pe13.nii.gz` 
  
  fslmaths fixation -Tmean -s 3 fixation_mean
  fslmaths effects -Tmean -s 3 effects_mean
  fslmaths derivatives -Tmean -s 3 derivatives_mean
  fslmaths stretch -Tmean -s 3 stretch_mean
  
  echo "Total volumes : " `fslnvols effects.nii.gz`
  rm -r run*/model*
   # Get sigma of residuals
#  fslmaths sigmasquared.nii.gz -sqrt sigma.nii.gz

# Prepare surfaces
if [ -e $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_ses-${SESSION}_T1w_brainmask.nii.gz ]; then
mri_vol2vol --mov $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_ses-${SESSION}_T1w_brainmask.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/effects.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask.nii.gz --nearest
else 
mri_vol2vol --mov $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/anat/sub-${SUBJECT}_T1w_brainmask.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/effects.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask.nii.gz --nearest
fi
fslmaths $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/effects.nii.gz -abs -Tmin -thr 0 -mas $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask -bin $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask

ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.inflated $SUBDIR/surf/lh.inflated
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.inflated $SUBDIR/surf/rh.inflated
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.curv $SUBDIR/surf/lh.curv
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.curv $SUBDIR/surf/rh.curv
ln -s $SUBJECTS_DIR/fsaverage/surf/lh.inflated $SUBDIR/surf/lh.fsaverage.inflated
ln -s $SUBJECTS_DIR/fsaverage/surf/rh.inflated $SUBDIR/surf/rh.fsaverage.inflated
ln -s $SUBJECTS_DIR/fsaverage/surf/lh.curv $SUBDIR/surf/lh.fsaverage.curv
ln -s $SUBJECTS_DIR/fsaverage/surf/rh.curv $SUBDIR/surf/rh.fsaverage.curv

# downsample surfaces
rm $SUBDIR/surf/*.ds.*


for hemi in rh lh; do
    for surface in pial white; do
        mri_surf2surf --hemi $hemi --srcsubject sub-${SUBJECT} --sval-xyz $surface --trgsubject fsaverage6 --trgicoorder 6 --trgsurfval  $SUBDIR/surf/${hemi}.${surface}.ds.mgh --tval-xyz $SUBJECTS_DIR/sub-${SUBJECT}/mri/T1.mgz 
        mris_convert --to-scanner $SUBDIR/surf/${hemi}.${surface}.ds.mgh $SUBDIR/surf/${hemi}.${surface}.ds.gii
done
        surface=inflated
        mri_surf2surf --hemi $hemi --srcsubject sub-${SUBJECT} --sval-xyz $surface --trgsubject fsaverage6 --trgicoorder 6 --trgsurfval  $SUBDIR/surf/${hemi}.${surface}.ds.mgh --tval-xyz $SUBJECTS_DIR/sub-${SUBJECT}/mri/T1.mgz
        mris_convert --to-scanner $SUBDIR/surf/${hemi}.${surface}.ds.mgh $SUBDIR/surf/${hemi}.${surface}.ds.gii
done
rm $SUBDIR/surf/*.ds.mgh


fi # PHASE 7

###############################################
# Run searchlight analysis
###############################################

if [ $PHASE == 0 ] || [ $PHASE == 8 ]; then

# Prepare engine
for hemi in rh lh; do
    python ~/Software/LeftHand/process/prepare_queryengine.py \
    $SUBDIR \
    $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/sequences.csv \
    $hemi $RADIUS "$RUNS"
done

fi #PHASE 8


TESTDIR=metrics
METRICS="acc_svm acc_svm_PCA invcompactness_correlation"
METRICS="acc_svm invcompactness_correlation"


if [ $PHASE == 0 ] || [ $PHASE == 9 ]; then
ACC_FWHM=5

# Run searchlight
for hemi in rh lh; do
    python ~/Software/LeftHand/process/surface_searchlight.py \
    $SUBDIR \
    $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/sequences.csv  \
    $SUBJECTS_DIR/sub-${SUBJECT}/label/${hemi}.cortex.label\
    $hemi $RADIUS $NPROC $TESTDIR "$RUNS"

# Resample to average
for metric in $METRICS; do
  mri_convert $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.gii $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.mgh
  ln -s $SUBJECTS_DIR/fsaverage6/surf/lh.sphere.reg $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.sphere6.reg
  ln -s $SUBJECTS_DIR/fsaverage6/surf/rh.sphere.reg $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.sphere6.reg 
  mri_surf2surf --srcsubject sub-${SUBJECT} \
    --srcsurfval $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.mgh \
    --trgsubject fsaverage \
    --srcsurfreg sphere6.reg \
    --trgsurfval $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.mgh \
    --hemi $hemi \
    --fwhm-trg $ACC_FWHM \
    --cortex

  mri_convert $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.mgh $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.gii
  done
done

rm $SUBDIR/surf/$TESTDIR/*.mgh

fi #PHASE 9

# Visualize

#freeview -f rh.inflated:overlay=rh.sl_accuracy_10.0.func.gii lh.inflated:overlay=lh.sl_accuracy_10.0.func.gii
#freeview -f rh.fsaverage.inflated:overlay=rh.sl_accuracy_10.0.fsaverage.func.gii:curv=rh.fsaverage.curv \
#lh.fsaverage.inflated:overlay=lh.sl_accuracy_10.0.fsaverage.func.gii:curv=lh.fsaverage.curv
do_show(){

OVERLAY_L=$1
OVERLAY_R=$2
THR_L=$3
THR_H=$4
FIG=$5

freeview -f $SURF_L:overlay=$OVERLAY_L:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
-viewport 3D  -colorscale -cam azimuth 180 elevation 0 -zoom 1.3 -ss ${FIG}_L_medial

freeview -f $SURF_L:overlay=$OVERLAY_L:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
-viewport 3D  -colorscale -cam azimuth 0 elevation 0 -zoom 1.2 -ss ${FIG}_L_lateral

freeview -f $SURF_R:overlay=$OVERLAY_R:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
-viewport 3D  -colorscale -cam azimuth 180 elevation 0 -zoom 1.2 -ss ${FIG}_R_lateral

freeview -f $SURF_R:overlay=$OVERLAY_R:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
-viewport 3D  -colorscale -cam azimuth 0 elevation 0 -zoom 1.3 -ss ${FIG}_R_medial

pngappend ${FIG}_L_medial.png + ${FIG}_R_medial.png ${FIG}_medial.png
pngappend ${FIG}_L_lateral.png + ${FIG}_R_lateral.png ${FIG}_lateral.png
pngappend ${FIG}_medial.png - ${FIG}_lateral.png ${FIG}.png
rm ${FIG}*lateral* ${FIG}*medial*

}

export SURF_L=$SUBJECTS_DIR/fsaverage/surf/lh.inflated
export SURF_R=$SUBJECTS_DIR/fsaverage/surf/rh.inflated

if [ $PHASE == 0 ] || [ $PHASE == 10 ]; then


for metric in $METRICS; do
      if [ $metric == 'invcompactness_correlation' ]; then
        LIM1=1
        LIM2=1.15
      else
        LIM1=0.30
        LIM2=0.40
      fi
      
      do_show $SUBDIR/surf/$TESTDIR/lh.sl_${metric}_${RADIUS}.fsaverage.func.gii $SUBDIR/surf/$TESTDIR/rh.sl_${metric}_${RADIUS}.fsaverage.func.gii $LIM1 $LIM2 $HOMEDIR/results/sl_${metric}_${RADIUS}_sub-${SUBJECT}_ses-${SESSION}

done

fi #PHASE 10

# average what we have
if [ $PHASE == 0 ] || [ $PHASE == 11 ]; then

 
for metric in $METRICS; do
    for hemi in rh lh; do
      maps=`echo $HOMEDIR/fmriprep/analysis/sub-*/ses-*/surf/metrics/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.gii`
      mris_calc -o $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii `echo $maps| cut -d" " -f1` mul 0
      echo $maps
      for mymap in $maps; do     
         mris_calc -o $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii add $mymap
      done
      mris_calc -o $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii div `echo $maps | wc -w`
    done

      if [ $metric == 'invcompactness_correlation' ]; then
#        mris_calc -o $HOMEDIR/fmriprep/analysis/surf/lh.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/surf/lh.sl_${metric}_${RADIUS}.mean.func.gii mul -1
#        mris_calc -o $HOMEDIR/fmriprep/analysis/surf/rh.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/surf/rh.sl_${metric}_${RADIUS}.mean.func.gii mul -1

        LIM1=1
        LIM2=1.15
      else
        LIM1=0.30
        LIM2=0.40
      fi
    
      do_show $HOMEDIR/fmriprep/analysis/surf/lh.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/surf/rh.sl_${metric}_${RADIUS}.mean.func.gii $LIM1 $LIM2 $HOMEDIR/results/sl_${metric}_${RADIUS}_mean

done

fi #PHASE 11

###############################################
# Extract roi data
###############################################



