#/bin/bash

# fix fieldmaps
# use mrqc instead of mriprep to get good segmentation
# try without using derivatives in the regressors...

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

# PHASE = 0 do everything

# put labels in $HOMEDIR/labels

# add $WD/freesurfer/fsaverage6
###############################################
# Variable definitions
###############################################
SUBJECT=$1
SESSION=$2
RESPONSES_FILE=$3 
export NTRIALS=$4 
PHASE=$5
RUNS=$6
RESPONSES_OPTIONS=$7
CONFOUND_INDICES=$8 
HOMEDIR=$9

OVERWRITE=0
ORGANIZE_FILE=select_organize_responses.py
PROGDIR=~/Software/LeftHand/process/
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/
#LABELSDIR=$HOMEDIR/labels

# Same always
SIMPLE3_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_3stretch.fsf
SIMPLE4_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_4stretch.fsf
SIMPLE6_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_6stretch_training.fsf
MULTI_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_multistretch.fsf

WD=$HOMEDIR/fmriprep
SUBJECTS_DIR=$WD/freesurfer
WORK=$HOMEDIR/work_sub-${SUBJECT}_ses-${SESSION}
ECHOTIME1=0.001544 #1551
ECHOTIME2=0.002544 
TE=25 #ms
EPIFAC=39
WFS=23.761 #2.366
SENSEP=3
FIELDSTRENGTH=7
EES=`echo "((1000 * $WFS)/($FIELDSTRENGTH*3.4*42.57 * ($EPIFAC+1))/$SENSEP)" | bc -l | awk '{printf "%f", $0}'`

FX_FILE=tstat
RADIUS=15.0
NPROC=30
MAXFEATPROCS=25
NEVS=4
ACC_FWHM=5
FWHM=5
SIGMA=`echo $FWHM/2.3548 | bc -l`

TESTDIR=metrics
METRICS="acc_svm acc_svm_PCA spread_correlation"
METRICS="within_spread_correlation spread_correlation acc_svm"
#METRICS="acc_svm spread_correlation within_spread_correlation"


if [ ! -e $WD ]; then 
    mkdir $WD
fi

# Find available fMRI runs
if [ "$RUNS" -eq "0" ]; then
    cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/
    RUNS0=`ls sub-${SUBJECT}_ses-${SESSION}_*bold.nii.gz  | cut -d'_' -f4 | cut -d'-' -f2 | cut -d'0' -f2 |tr '\r\n' ' '| sed '$s/ $/\n/g'`
    RUNS=""
      for run in $RUNS0; do
         if [ -e $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run ]; then
           RUNS="$RUNS $run"
         fi
      done  
    RUNS=`echo $RUNS | sed 's/^[\t ]*//g'`
    echo "Found following runs: $RUNS"
fi

###############################################
# Create explanatory variables for fMRI analysis
###############################################
if [ $PHASE == 0 ] || [ $PHASE == 1 ]; then

mkdir -p $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}
cd $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}
rm -r $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run*
python $PROGDIR/$ORGANIZE_FILE $RESPONSES_FILE ${SUBJECT} ${SESSION} $RESPONSES_OPTIONS > CorrectTrials.csv

fi


###############################################
# Prepare fieldmaps and structural images
###############################################
if [ $PHASE == 0 ] || [ $PHASE == 2 ]; then

cd $HOMEDIR/

cp data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_epi1.nii.gz data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_magnitude.nii.gz 

fslmaths data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_epi2.nii.gz -mul 6.28 data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap_rads.nii.gz 
fugue --loadfmap=data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap_rads.nii.gz -m --savefmap=data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap.nii.gz

# fix Philips json files
for run in $RUNS; do
    sed -i 's/Axis/Direction/g' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json
    myline=`cat data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json | grep EchoTrainLength`
    sed -i '/'"$myline"'/a \  \"EffectiveEchoSpacing\": '"$EES"',' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json    
done 

cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/
INTENDEDFOR=\"`echo ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-??_bold.nii.gz | sed 's/ /","/g'`\"
echo -e "{\n\"EchoTime1\": $ECHOTIME1, \n\"EchoTime2\": $ECHOTIME2, \n\"PhaseEncodingDirection\": \"j\", \n\"Units\": \"rad/s\", \n\"IntendedFor\": [ $INTENDEDFOR ] \n}" > $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap.json 

cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/
ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec2.nii.gz re.nii.gz 
ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.nii.gz im.nii.gz

fslroi sub-${SUBJECT}_ses-${SESSION}_MP2RAGE.nii.gz mag 1 1
#fslroi sub-${SUBJECT}_ses-${SESSION}_MP2RAGE.nii.gz phase 0 1

cd $PROGDIR/mprageconvert/
matlab -nosplash -nodisplay -r "addpath $PROGDIR/mprageconvert/Nifti_tools; create_mp2rage_command('$HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/', 're.nii.gz ', 'im.nii.gz', 'sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.json'); exit"

cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/
fslmaths MP2RAGE.nii.gz -add 0.5 MP2RAGEpos.nii.gz 
fslmaths mag -mul MP2RAGEpos magRAGE
 
fi
###############################################
# Create singularity file container
###############################################

if [ $PHASE == 0 ] || [ $PHASE == 3 ]; then

if [ ! -e $HOMEDIR/fmriprep.img ]; then
    cd /tmp
    docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /tmp/fmriprep1.4:/output \
        singularityware/docker2singularity \
        poldracklab/fmriprep:1.4.0
        
#    docker run --privileged -t --rm \
#        -v /var/run/docker.sock:/var/run/docker.sock \
#        -v /home/benjamin.garzon/Data/LeftHand/Lund1test/fmriprep1.3:/output \
#        singularityware/docker2singularity \
#        poldracklab/fmriprep:1.3.0.post1
    
    cp /tmp/fmriprep/poldracklab_fmriprep_latest-*.img $HOMEDIR/fmriprep.img
fi



###############################################
# Process with fmriprep
###############################################


rm -r $WORK
mkdir $WORK
cp /usr/local/freesurfer/license.txt $WD/

if [ ]; then

PYTHONPATH="" singularity run \
   $HOMEDIR/fmriprep.img \
   $HOMEDIR/data_BIDS \
   $WD \
   participant \
   --ignore slicetiming \
   --ignore fieldmaps \
   --longitudinal \
   --fs-license-file \
   $WD/license.txt \
   --nthreads 4 \
   -w $WORK \
   --write-graph \
   --output-space 'fsnative' \
   --participant-label ${SUBJECT}   

fi

PYTHONPATH="" singularity run \
   $HOMEDIR/fmriprep.img \
   $HOMEDIR/data_BIDS \
   $WD \
   participant \
   --ignore slicetiming \
   --ignore fieldmaps \
   --longitudinal \
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
   --longitudinal \
   --fs-license-file \
   $WD/license.txt \
   --nthreads 4 \
   -w $WORK \
   --write-graph \
   --template 'MNI152NLin2009cAsym' \
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
mkdir -p $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/
mkdir $SUBDIR
mkdir $SUBDIR/surf

# collect the data
fslmerge -t $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data.nii.gz $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-*_bold_space-T1w_preproc.nii.gz
fslmaths $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data.nii.gz -Tmean $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data_mean.nii.gz
rm $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data.nii.gz
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

if [ "$OVERWRITE" -eq 1 ]; then
    for run in $RUNS; 
    do    
          RUN_DIR=$SUBDIR/run$run
          rm -r $RUN_DIR
    done
fi

# Compute activation maps
if [ $PHASE == 0 ] || [ $PHASE == 5 ]; then

cd $SUBDIR

for xxx in 1; do # just to make it wait
for run in $RUNS; 
do
  if [ ! -e "$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/run$run/volume" ]; then
      FUNCVOL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_space-T1w_preproc.nii.gz
      FUNCSURFL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_space-fsnative.L.func.gii
      FUNCSURFR=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_space-fsnative.R.func.gii
    
      RUN_DIR=$SUBDIR/run$run
    
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
    
      cp $SIMPLE6_FSF_FILE $FSF_FILE
      sed -i "s%@EV1%$RUN_DIR/EV1.csv%" $FSF_FILE        
      sed -i "s%@EV2%$RUN_DIR/EV2.csv%" $FSF_FILE        
      sed -i "s%@EV3%$RUN_DIR/EV3.csv%" $FSF_FILE        
      sed -i "s%@EV4%$RUN_DIR/EV4.csv%" $FSF_FILE   
      sed -i "s%@EV5%$RUN_DIR/EV5.csv%" $FSF_FILE        
      sed -i "s%@EV6%$RUN_DIR/EV6.csv%" $FSF_FILE        

      # if there are empty regressors, set them to 0
      for k in 3 4 5 6; do          
          if [ `cat EV${k}.csv | wc -l` == 0 ]; then
           sed -i "s%set fmri(convolve${k}) 3%set fmri(convolve${k}) 0%" $FSF_FILE 
           sed -i "s%set fmri(shape${k}) 3%set fmri(shape${k}) 10%" $FSF_FILE 
          fi
      done
      
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
  fi # volume exists
done

done
rm $SUBDIR/run*/smoothed_data.nii.gz
fi # PHASE 5

  ##### do second level analysis


# Prepare parameter estimates. Using fmriprep outputs
# assume PHASE5 already run
if [ $PHASE == 0 ] || [ $PHASE == 6 ]; then

if [ ! -e "$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/effects.nii.gz" ]; then
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
  feat $FSF_FILE &
  sleep 10
    while [ `ps -a | grep fsl_sub | wc -l` -gt $MAXFEATPROCS ]; do
      sleep 1800
    done 
  done 


# wait for all to finish
while [ `ls -v $RUN_DIR/model*/analysis.feat/stats/${FX_FILE}2.nii.gz | wc -l` -lt $NTRIALS ]; do
    echo "Waiting for all models to finish"
    echo "sub-${SUBJECT}: Only `ls -v $RUN_DIR/model*/analysis.feat/stats/${FX_FILE}2.nii.gz | wc -w` out of $NTRIALS have finished"
    sleep 1800
done

  fslmerge -t $RUN_DIR/effects_$run `ls -v $RUN_DIR/model*/analysis.feat/stats/${FX_FILE}1.nii.gz` 
  fslmerge -t $RUN_DIR/derivatives_$run `ls -v $RUN_DIR/model*/analysis.feat/stats/${FX_FILE}2.nii.gz` 
  rm -r $RUN_DIR/model*

done


#echo done
#fi # PHASE 6
#if [ $PHASE == 0 ] || [ $PHASE == 7 ]; then



# Merge parameters
cd $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/

#for run in $RUNS; do
#  fslmerge -t run$run/effects_$run `ls -v run$run/model*/analysis.feat/stats/${FX_FILE}1.nii.gz` 
#  fslmerge -t run$run/derivatives_$run `ls -v run$run/model*/analysis.feat/stats/${FX_FILE}2.nii.gz` 
#done
  fslmerge -t effects `ls -v run*/effects*` 
  fslmerge -t derivatives `ls -v run*/derivatives*` 

#  fslmerge -t other1 `ls -v run$run/model*/analysis.feat/stats/pe3.nii.gz` 
#  fslmerge -t other2fsaverage `ls -v run$run/model*/analysis.feat/stats/pe5.nii.gz` 
#  fslmerge -t other3 `ls -v run$run/model*/analysis.feat/stats/pe7.nii.gz` 
#  fslmerge -t other4 `ls -v run$run/model*/analysis.feat/stats/pe9.nii.gz` 

#  fslmerge -t fixation `ls -v run$run/model*/analysis.feat/stats/pe11.nii.gz` 
#  fslmerge -t stretch `ls -v run$run/model*/analysis.feat/stats/pe13.nii.gz` 
  
#  fslmaths fixation -Tmean -s 3 fixation_mean
  fslmaths effects -Tmean -s 3 effects_mean
  fslmaths derivatives -Tmean -s 3 derivatives_mean
#  fslmaths stretch -Tmean -s 3 stretch_mean
  
  echo "Total volumes : " `fslnvols effects.nii.gz`
  rm -r run*/model*
  rm fixation.nii.gz stretch.nii.gz run*/effects* run*/derivatives*
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

fi # check if effects already exists
fi # PHASE 7

###############################################
# Run searchlight analysis
###############################################

if [ $PHASE == 0 ] || [ $PHASE == 8 ]; then

# Prepare engine
for hemi in rh lh; do
    python $PROGDIR/prepare_queryengine.py \
    $SUBDIR \
    $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/sequences.csv \
    $hemi $RADIUS "$RUNS"
done

fi #PHASE 8


if [ $PHASE == 0 ] || [ $PHASE == 9 ]; then

# Run searchlight
for hemi in rh lh; do
    python $PROGDIR/surface_searchlight.py \
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

rm $SUBDIR/surf/$TESTDIR/*.mgh $SUBDIR/surf/$TESTDIR/qe.mask.nii.gz #$SUBDIR/surf/*gzipped.hdf5 

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

if [ ! -e $HOMEDIR/fmriprep/results ]; then 
    mkdir $HOMEDIR/fmriprep/results
fi

for metric in $METRICS; do
      if [ $metric == 'spread_correlation' ]  || [ $metric == 'within_spread_correlation' ]; then
        LIM1=1
        LIM2=1.20
      else
        LIM1=0.30
        LIM2=0.40
      fi
      
      do_show $SUBDIR/surf/$TESTDIR/lh.sl_${metric}_${RADIUS}.fsaverage.func.gii $SUBDIR/surf/$TESTDIR/rh.sl_${metric}_${RADIUS}.fsaverage.func.gii $LIM1 $LIM2 $HOMEDIR/fmriprep/results/sl_${metric}_${RADIUS}_sub-${SUBJECT}_ses-${SESSION}

done

fi #PHASE 10

# average what we have and find rois
if [ $PHASE == 0 ] || [ $PHASE == 11 ]; then

LABELSDIR=$HOMEDIR/fmriprep/analysis/sub-$SUBJECT/labels

if [ ! -e $HOMEDIR/fmriprep/analysis/surf ]; then 
    mkdir $HOMEDIR/fmriprep/analysis/surf
fi
if [ ! -e $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf ]; then 
    mkdir $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf
    mkdir $LABELSDIR
fi


# get labels
#labels of interest
LABELS="S_front_inf S_postcentral G_precentral S_front_sup Pole_occipital S_intrapariet_and_P_trans G_temporal_middle G_and_S_cingul-Mid-Ant"

cd $LABELSDIR

mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi rh --annotation aparc.a2009s --subject fsaverage6
mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi lh --annotation aparc.a2009s --subject fsaverage6

for label in $LABELS; do
    echo $label
    mris_convert --label lh.${label}.label ${label}\
        $SUBJECTS_DIR/fsaverage6/surf/lh.white \
        lh.${label}.label.gii
    $HCPDIR/wb_command -gifti-label-to-roi lh.${label}.label.gii -name ${label} lh.${label}.func.gii
        
    mris_convert --label rh.${label}.label ${label}\
        $SUBJECTS_DIR/fsaverage6/surf/rh.white \
        rh.${label}.label.gii
    $HCPDIR/wb_command -gifti-label-to-roi rh.${label}.label.gii -name ${label} rh.${label}.func.gii

done

# done every session, although it would be enough to do it for the last one
for metric in $METRICS; do
# average across subjects
    for hemi in rh lh; do
      maps=`echo $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/ses-*/surf/metrics/$hemi.sl_${metric}_${RADIUS}.func.gii`
      echo $maps
      rm $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii

      $HCPDIR/wb_command -metric-math "(x*0)" $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii -var x `echo $maps| cut -d" " -f1` -column 2
      for mymap in $maps; do  
      $HCPDIR/wb_command -metric-math "(x+y)" $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii -var x $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii -var y $mymap -column 2

      done
      mris_calc -o $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii div `echo $maps | wc -w`
    done #hemi

    
      if [ $metric == 'spread_correlation' ] || [ $metric == 'within_spread_correlation' ]; then

        LIM1=1
        LIM2=1.20
      else
        LIM1=0.30
        LIM2=0.40
      fi

      cd $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/
      echo HEMI NAME INDEXMAX > $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/lh.${metric}_maxima.csv
      echo HEMI NAME INDEXMAX > $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/rh.${metric}_maxima.csv
      
      for label in $LABELS; do
      
      $HCPDIR/wb_command -metric-math "(x*y)" $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/aux.func.gii \
      -var x $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/lh.sl_${metric}_${RADIUS}.mean.func.gii \
      -var y $LABELSDIR/lh.${label}.func.gii      
      echo lh $label `$HCPDIR/wb_command -metric-stats $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/aux.func.gii -reduce INDEXMAX` >> $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/lh.${metric}_maxima.csv

      $HCPDIR/wb_command -metric-math "(x*y)" $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/aux.func.gii \
      -var x $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/rh.sl_${metric}_${RADIUS}.mean.func.gii \
      -var y $LABELSDIR/rh.${label}.func.gii

      echo rh $label `$HCPDIR/wb_command -metric-stats $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/aux.func.gii -reduce INDEXMAX` >> $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/rh.${metric}_maxima.csv
      done

      rm $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/aux.func.gii

# generate results in fsaverage space
#average across sessions
    for hemi in rh lh; do
      maps=`echo $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/ses-*/surf/metrics/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.gii`
      mris_calc -o $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii `echo $maps| cut -d" " -f1` mul 0
      echo $maps
      for mymap in $maps; do     
         mris_calc -o $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii add $mymap
      done #maps
      mris_calc -o $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/$hemi.sl_${metric}_${RADIUS}.mean.func.gii div `echo $maps | wc -w`
    done #hemi
      
      do_show $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/lh.sl_${metric}_${RADIUS}.mean.func.gii \
      $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/rh.sl_${metric}_${RADIUS}.mean.func.gii \
      $LIM1 $LIM2 $HOMEDIR/fmriprep/results/sl_${metric}_${RADIUS}_sub-${SUBJECT}_mean
      
#average across subjects      
    for hemi in rh lh; do
      maps=`echo $HOMEDIR/fmriprep/analysis/sub-*/ses-${SESSION}/surf/metrics/$hemi.sl_${metric}_${RADIUS}.fsaverage.func.gii`
      mris_calc -o $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.ses-${SESSION}.mean.func.gii `echo $maps| cut -d" " -f1` mul 0
      echo $maps
      for mymap in $maps; do     
         mris_calc -o $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.ses-${SESSION}.mean.func.gii $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.ses-${SESSION}.mean.func.gii add $mymap
      done #maps
      mris_calc -o $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.ses-${SESSION}.mean.func.gii $HOMEDIR/fmriprep/analysis/surf/$hemi.sl_${metric}_${RADIUS}.ses-${SESSION}.mean.func.gii div `echo $maps | wc -w`
    done #hemi
      
      do_show $HOMEDIR/fmriprep/analysis/surf/lh.sl_${metric}_${RADIUS}.ses-${SESSION}.mean.func.gii \
      $HOMEDIR/fmriprep/analysis/surf/rh.sl_${metric}_${RADIUS}.ses-${SESSION}.mean.func.gii \
      $LIM1 $LIM2 $HOMEDIR/fmriprep/results/sl_${metric}_${RADIUS}_ses-${SESSION}_mean

done

rm $LABELSDIR/*label*
fi #PHASE 11

###############################################
# Extract roi data
###############################################
# average what we have
if [ $PHASE == 0 ] || [ $PHASE == 12 ]; then
metric=spread_correlation
mkdir $SUBDIR/surf/rois
for hemi in rh lh; do
    python $PROGDIR/spread_from_rois.py \
    $SUBDIR \
    $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/sequences.csv \
    $SUBJECT \
    $SESSION \
    $hemi \
    $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf/${hemi}.${metric}_maxima.csv \
    $RADIUS $NPROC $TESTDIR "$RUNS" correlation   
done

    
# gather all the data together
cat $HOMEDIR/fmriprep/analysis/sub-*/ses-*/surf/metrics/*roi_correlation_distances.csv > $HOMEDIR/fmriprep/analysis/surf/roi_correlation_distances.csv
fi 


