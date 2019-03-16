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

###############################################
# Variable definitions
###############################################
SUBJECT=105 #$1
SESSION=1
export NTRIALS=28 # $2
RESPONSES_FILE=trialsfile-lup2s007_fmri.csv

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lundpilot1
RESPONSES_OPTIONS="2.8 0.5 0 25 4"
PROGDIR=~/Software/LeftHand/process/
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

# Same always
SIMPLE_FSF_FILE=$HOMEDIR/fMRInoreg_4stretch.fsf
MULTI_FSF_FILE=$HOMEDIR/fMRInoreg_multistretch.fsf
RUNS="1 2 3 4"
WD=$HOMEDIR/fmriprep
SUBJECTS_DIR=$WD/freesurfer
WORK=$HOMEDIR/work${SUBJECT}
ECHOTIME1=0.001544
ECHOTIME2=0.002544 
TE=25 #ms
EPIFAC=39
WFS=22.366
SENSEP=3
FIELDSTRENGTH=7
EES=`echo "((1000 * $WFS)/($FIELDSTRENGTH*3.4*42.57 * ($EPIFAC+1))/$SENSEP)" | bc -l | awk '{printf "%f", $0}'`

INTENDEDFOR="[\"func/sub-${SUBJECT}_task-sequence_run-01_bold.nii.gz\", \"func/sub-${SUBJECT}_task-sequence_run-02_bold.nii.gz\", \"func/sub-${SUBJECT}_task-sequence_run-03_bold.nii.gz\", \"func/sub-${SUBJECT}_task-sequence_run-04_bold.nii.gz\"]"
CONFOUND_INDICES=28,29,30,31,32,33
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
cd $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}
rm -r run*
python $PROGDIR/organize_responses.py $RESPONSES_FILE $RESPONSES_OPTIONS


###############################################
# Prepare fieldmaps
###############################################
cd $HOMEDIR/
mv data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_epi1.nii.gz data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_magnitude.nii.gz 

fslmaths data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_epi2.nii.gz -mul 6.28 data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_fieldmap_rads.nii.gz 
fugue --loadfmap=data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_fieldmap_rads.nii.gz -m --savefmap=data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_fieldmap.nii.gz

for run in $RUNS; do
    sed -i 's/Axis/Direction/g' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold.json
    myline=`cat data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold.json | grep EchoTrainLength`
    sed -i '/'"$myline"'/a \  \"EffectiveEchoSpacing\": '"$EES"',' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold.json    
done 

echo -e "{\n\"EchoTime1\": $ECHOTIME1, \n\"EchoTime2\": $ECHOTIME2, \n\"PhaseEncodingDirection\": \"j\", \n\"Units\": \"rad/s\", \n\"IntendedFor\": $INTENDEDFOR\n}" > data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_fieldmap.json 

###############################################
# Create singularity file container
###############################################

if [ ! -e $HOMEDIR/fmriprep.img ]; then
    cd /tmp
    docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /tmp/fmriprep:/output \
        singularityware/docker2singularity \
        poldracklab/fmriprep:latest
    
    cp /tmp/fmriprep/poldracklab_fmriprep_latest-*.img $HOMEDIR/fmriprep.img
fi

rm -r $WORK
mkdir $WORK
cp /usr/local/freesurfer/license.txt $WD/


###############################################
# Process with fmriprep
###############################################
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
   
# --ignore fieldmaps \
#   --use-syn-sdc \
#   --use-aroma \

# collect the data
fslmerge -t $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data.nii.gz $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_task-sequence_run-*_bold_space-T1w_preproc.nii.gz
fslmaths $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data.nii.gz -Tmean $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data_mean.nii.gz


###############################################
# Run fMRI analysis
###############################################
# Native surface space from fmriprep
#STRUCT=$HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_T1w_brain.nii.gz

SUBDIR=$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/
mkdir $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/
mkdir $SUBDIR
mkdir $SUBDIR/surf
cd $SUBDIR

# To facilitate visualization
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.pial $SUBDIR/surf/lh.pial.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.pial $SUBDIR/surf/rh.pial.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.white $SUBDIR/surf/lh.white.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.white $SUBDIR/surf/rh.white.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.midthickness $SUBDIR/surf/lh.midthickness.gii
mris_convert --to-scanner $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.midthickness $SUBDIR/surf/rh.midthickness.gii
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.sphere.reg $SUBDIR/surf/lh.sphere.reg
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.sphere.reg $SUBDIR/surf/rh.sphere.reg

for run in $RUNS; 
do

  FUNCVOL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold_space-T1w_preproc.nii.gz
  FUNCSURFL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold_space-fsnative.L.func.gii
  FUNCSURFR=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold_space-fsnative.R.func.gii

  #Smoothing 
  fslmaths $FUNCVOL -s $FWHM $RUN_DIR/smoothed_data.nii.gz
  mri_surf2surf --s sub-${SUBJECT} --sval $FUNCSURFL --hemi lh --fwhm-trg $FWHM --tval $RUN_DIR/lh.smoothed.func.gii
  mri_surf2surf --s sub-${SUBJECT} --sval $FUNCSURFR --hemi rh --fwhm-trg $FWHM --tval $RUN_DIR/rh.smoothed.func.gii
  $HCPDIR/wb_command -metric-dilate $RUN_DIR/lh.smoothed.func.gii $SUBDIR/surf/lh.midthickness.gii 100 $RUN_DIR/lh.smoothed.func.gii -nearest 
  $HCPDIR/wb_command -metric-dilate $RUN_DIR/rh.smoothed.func.gii $SUBDIR/surf/rh.midthickness.gii 100 $RUN_DIR/rh.smoothed.func.gii -nearest 

  RUN_DIR=$SUBDIR/run$run
#  rm -r $RUN_DIR
#  mkdir $RUN_DIR
  cut -d$'\t' -f $CONFOUND_INDICES $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold_confounds.tsv | tail -n +2 > $RUN_DIR/mc_run-0${run}.csv

  FSF_FILE=$RUN_DIR/fMRI.fsf
  
  cp $SIMPLE_FSF_FILE $FSF_FILE
  cp $HOMEDIR/responses/sub-${SUBJECT}/run$run/*.csv $MODEL_DIR/$RUN_DIR/
  
  cd $RUN_DIR

  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@STRUCTURAL_BRAIN%$STRUCT%" $FSF_FILE
  sed -i "s%@ANALYSIS%$RUN_DIR/analysis%" $FSF_FILE  
  sed -i "s%@CONFOUNDS%$RUN_DIR/mc_run-0${run}.csv%" $FSF_FILE        
  
  for EV in `seq $NEVS`; 
  do  
     sed -i "s%@EV$EV%$RUN_DIR/EV$EV.csv%" $FSF_FILE        
  done

  # now run it
  feat_model fMRI #$RUN_DIR/mc_run-0${run}.csv
  rm -r $RUN_DIR/surfL $RUN_DIR/surfR $RUN_DIR/volume 
#  film_gls --rn=$RUN_DIR/surfL --sa --in=$RUN_DIR/lh.smoothed.func.gii --pd=fMRI.mat --con=fMRI.con --mode=surface --in2=$SUBDIR/surf/lh.midthickness.gii
#  film_gls --rn=$RUN_DIR/surfR --sa --in=$RUN_DIR/rh.smoothed.func.gii --pd=fMRI.mat --con=fMRI.con --mode=surface --in2=$SUBDIR/surf/rh.midthickness.gii 
  film_gls --rn=$RUN_DIR/volume --sa --ms=$FWHM --in=$RUN_DIR/smoothed_data.nii.gz --thr=1000 --pd=fMRI.mat --con=fMRI.con --mode=volume 

done


# Using fmriprep outputs

for run in $RUNS; 
do

  FUNC=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/func/sub-${SUBJECT}_task-sequence_run-0${run}_bold_space-T1w_preproc.nii.gz

  RUN_DIR=$SUBDIR/run$run
#  rm -r $RUN_DIR
  mkdir $RUN_DIR 

    # Get confounds 
  
  for TRIAL in `seq $NTRIALS`; 
    do
   echo $TRIAL

  MODEL_DIR=$RUN_DIR/model$TRIAL
  mkdir $MODEL_DIR
  FSF_FILE=$MODEL_DIR/fMRI.fsf
  
  cp $MULTI_FSF_FILE $FSF_FILE
  cp $HOMEDIR/responses/sub-${SUBJECT}/run$run/TRIAL$TRIAL $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/run$run/OTHER$TRIAL $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/run$run/FIXATION $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/run$run/STRETCH $MODEL_DIR/
  
  cd $MODEL_DIR

  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@ANALYSIS%$MODEL_DIR/analysis%" $FSF_FILE
  
  sed -i "s%@EV1%$MODEL_DIR/TRIAL$TRIAL%" $FSF_FILE        
  sed -i "s%@EV2%$MODEL_DIR/OTHER$TRIAL%" $FSF_FILE        
  sed -i "s%@EV3%$MODEL_DIR/FIXATION%" $FSF_FILE        
  sed -i "s%@EV4%$MODEL_DIR/STRETCH%" $FSF_FILE        
  sed -i "s%@CONFOUNDS%$RUN_DIR/mc_run-0${run}.csv%" $FSF_FILE        

  # Now run it
  feat $FSF_FILE 
  done &
done

# Merge parameters
cd $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}
for run in $RUNS; do
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
  
# Get sigma of residuals
fslmaths sigmasquared.nii.gz -sqrt sigma.nii.gz


###############################################
# Run searchlight analysis
###############################################

# Prepare surfaces
mri_vol2vol --mov $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/anat/sub-${SUBJECT}_T1w_brainmask.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/effects.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/mask.nii.gz --nearest
fslmaths $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/effects.nii.gz -abs -Tmin -thr 0 -mas $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/mask -bin $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/mask

mkdir $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf

ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.inflated $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/lh.inflated
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.inflated $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/rh.inflated
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/lh.curv $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/lh.curv
ln -s $SUBJECTS_DIR/sub-${SUBJECT}/surf/rh.curv $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/rh.curv
ln -s $SUBJECTS_DIR/fsaverage/surf/lh.inflated $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/lh.fsaverage.inflated
ln -s $SUBJECTS_DIR/fsaverage/surf/rh.inflated $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/rh.fsaverage.inflated
ln -s $SUBJECTS_DIR/fsaverage/surf/lh.curv $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/lh.fsaverage.curv
ln -s $SUBJECTS_DIR/fsaverage/surf/rh.curv $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/surf/rh.fsaverage.curv


# Prepare engine
for hemi in rh lh; do
    python ~/Software/LeftHand/process/prepare_queryengine.py \
    $HOMEDIR/fmriprep/analysis/sub-${SUBJECT} \
    $HOMEDIR/responses/sub-${SUBJECT}/sequences.csv \
    $hemi $RADIUS $NPROC
done

TESTDIR=with_PCA
METRICS="acc_svm correlation"
ACC_FWHM=5

# Run searchlight
for hemi in rh lh; do
    python ~/Software/LeftHand/process/surface_searchlight.py \
    $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}  \
    $HOMEDIR/responses/sub-${SUBJECT}/sequences.csv \
    $SUBJECTS_DIR/sub-${SUBJECT}/label/${hemi}.cortex.label\
    $hemi $RADIUS $NPROC

# Resample to average
for metric in $METRICS; do
  mri_convert $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.gii $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.mgh

  mri_surf2surf --srcsubject sub-${SUBJECT} \
    --srcsurfval $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.mgh \
    --trgsubject fsaverage \
    --trgsurfval $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.mgh \
    --hemi $hemi \
    --fwhm-trg $ACC_FWHM \
    --cortex

  mri_convert $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.mgh $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.gii
  done
done

rm $SUBDIR/surf/$TESTDIR/*.mgh

# Visualize

#freeview -f rh.inflated:overlay=rh.sl_accuracy_10.0.func.gii lh.inflated:overlay=lh.sl_accuracy_10.0.func.gii
#freeview -f rh.fsaverage.inflated:overlay=rh.sl_accuracy_10.0.fsaverage.func.gii:curv=rh.fsaverage.curv \
#lh.fsaverage.inflated:overlay=lh.sl_accuracy_10.0.fsaverage.func.gii:curv=lh.fsaverage.curv

mkdir $SUBDIR/results


do_show(){

OVERLAY_L=$1
OVERLAY_R=$2
THR_L=$3
THR_H=$4
FIG=$5

#freeview -f $SURF_L:overlay=$OVERLAY_L:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
#-viewport 3D -colorscale -cam azimuth 180 elevation 0 -zoom 1.3 -ss ${FIG}_L_medial

#freeview -f $SURF_L:overlay=$OVERLAY_L:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
#-viewport 3D  -colorscale -cam azimuth 0 elevation 0 -zoom 1 -ss ${FIG}_L_lateral

freeview -f $SURF_R:overlay=$OVERLAY_R:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
-viewport 3D  -colorscale -cam azimuth 180 elevation 0 -zoom 1.3 -ss ${FIG}_R_lateral

freeview -f $SURF_R:overlay=$OVERLAY_R:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse' \
-viewport 3D  -colorscale -cam azimuth 0 elevation 0 -zoom 1 -ss ${FIG}_R_medial

}

export SURF_L=$SUBJECTS_DIR/fsaverage/surf/lh.inflated
export SURF_R=$SUBJECTS_DIR/fsaverage/surf/rh.inflated

do_show $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.gii $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.gii 0 1 $SUBDIR/results/sl_${metric}_${RADIUS}


###############################################
# Extract roi data
###############################################



