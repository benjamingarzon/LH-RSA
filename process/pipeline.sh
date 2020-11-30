#/bin/bash

# PHASE 0 - run everything
# PHASE 1 - create EV
# PHASE 2 - prepare fieldmaps and structural images
# PHASE 3 - Process with fMRIprep 
# PHASE 4 - Prepare surfaces, labels and masks
# PHASE 5 - Run fMRI analysis
# PHASE 6 - Obtain single-trial parameter estimates
# PHASE 7 - Prepare engine for searchlight analysis
# PHASE 8 - Searchlight analysis
# PHASE 9 - Plot surface results


###############################################################################
# Variable definitions
###############################################################################
# parameters
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
HCPDIR=~/Software/workbench/bin_rh_linux64/
PATH=$PATH:$HCPDIR
INVALID_FILES=$HOMEDIR/logs/invalid_files.txt
WD=$HOMEDIR/fmriprep
SUBJECTS_DIR=$WD/freesurfer
WORK=$HOMEDIR/work_sub-${SUBJECT}_ses-${SESSION}

# fMRI design files
SIMPLE6_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_6stretch_training_noF.fsf
MULTI_FSF_FILE=$HOMEDIR/fmri_designs/fMRInoreg_multistretch.fsf
SINGLE1_FSF_FILE=$HOMEDIR/fmri_designs/fMRIsingle1.fsf
SINGLE2_FSF_FILE=$HOMEDIR/fmri_designs/fMRIsingle2.fsf

# MR parameters
ECHOTIME1=0.001544 
ECHOTIME2=0.002544 
TE=25 #ms
EPIFAC=39
WFS=23.761 
SENSEP=3
FIELDSTRENGTH=7
TR=1.2
TemporalFilter=90
RI=-500
RS=0.2442

# Analysis parameters
FX_FILE=pe #tstat
RADIUS=10.0
NEVS=4
ACC_FWHM=10
FWHM=8
FWHM_SURF=10
SIGMA=`echo $FWHM/2.3548 | bc -l`
TESTDIR=metrics
METRICS="within_spread_correlation secmom acc_svm"

if [ ! -e $WD ]; then 
    mkdir $WD
fi

# Find available fMRI runs
if [ "$RUNS" -eq "0" ]; then
    cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/
    RUNS0=`ls sub-${SUBJECT}_ses-${SESSION}_*bold.nii.gz  | cut -d'_' -f4 | cut -d'-' -f2 | cut -d'0' -f2 |tr '\r\n' ' '| sed '$s/ $/\n/g'`
    RUNS=""
	# Select those with responses only
      for run in $RUNS0; do
        #if [ -e $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run ]; then
           RUNS="$RUNS $run"
        # fi
      done  
    RUNS=`echo $RUNS | sed 's/^[\t ]*//g'`
    echo "Found following runs: $RUNS"
fi

###############################################################################
# Create explanatory variables for fMRI analysis
###############################################################################
if [ $PHASE == 0 ] || [ $PHASE == 1 ]; then

mkdir -p $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}
cd $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}
rm -r $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run*
python $PROGDIR/$ORGANIZE_FILE $RESPONSES_FILE ${SUBJECT} ${SESSION} $RESPONSES_OPTIONS > CorrectTrials.csv

fi


###############################################################################
# Prepare fieldmaps and structural images
###############################################################################
if [ $PHASE == 0 ] || [ $PHASE == 2 ]; then

cd $HOMEDIR/

# fix Philips json files
for run in $RUNS; do
    sed -i 's/Axis/Direction/g' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json
    myline=`cat data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json | grep EchoTrainLength`
    sed -i '/EffectiveEchoSpacing/d' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json
    sed -i '/WaterFatShift/d' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json
    sed -i '/ParallelReductionFactorInPlane/d' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json
    sed -i '/'"$myline"'/a \  \"ParallelReductionFactorInPlane\": '"$SENSEP"',' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json    
    sed -i '/'"$myline"'/a \  \"WaterFatShift\": '"$WFS"',' data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json    

    #check that the runs are ok
    nvols=`fslnvols data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.nii.gz`
    echo "Run $run: $nvols vols"
    if [ "$nvols" -eq 1 ]; then
        echo "Found invalid run: data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.nii.gz"

	mv data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.nii.gz data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_invalid.nii.gz
	mv data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.json data_BIDS/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold_invalid.json
    fi

done 

# gather all invalid files
find data_BIDS | grep invalid.nii.gz > $INVALID_FILES

# create fieldmaps
if [ `ls data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_run*_epi1.nii.gz | wc -l` -gt 0 ]; then
    # more than one?
    for run in $RUNS; do
        cp data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_run-0${run}_epi1.nii.gz \
data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_run-0${run}_magnitude.nii.gz 
        fslmaths data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_run-0${run}_epi2.nii.gz \
-mul $RS -add $RI data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_run-0${run}_fieldmap.nii.gz 

        INTENDEDFOR="ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_bold.nii.gz"
        echo -e "{\n\"PhaseEncodingDirection\": \"j\", \n\"Units\": \"Hz\", \n\"IntendedFor\": [ \"$INTENDEDFOR\" ] \n}" > $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_run-0${run}_fieldmap.json 
    done
else
    # only one fieldmap
    cp data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_epi1.nii.gz data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_magnitude.nii.gz 

    # rescale values see PMC3998686
    # FP = (PV + RI/RS)/SS = PV*RS + RI
    fslmaths data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_epi2.nii.gz -mul $RS -add $RI data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap.nii.gz 

    cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/
    INTENDEDFOR=\"`echo ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-??_bold.nii.gz | sed 's/ /","/g'`\"
    echo -e "{\n\"PhaseEncodingDirection\": \"j\", \n\"Units\": \"Hz\", \n\"IntendedFor\": [ $INTENDEDFOR ] \n}" > $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/fmap/sub-${SUBJECT}_ses-${SESSION}_fieldmap.json 

fi

# structural data
cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/
if [ -e sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.nii.gz ]; then
    ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec2.nii.gz re.nii.gz 
    ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.nii.gz im.nii.gz
else
    ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRecReal.nii.gz re.nii.gz 
    ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRecImag.nii.gz im.nii.gz
    ln -s sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRecImag.json sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.json
fi

fslroi sub-${SUBJECT}_ses-${SESSION}_MP2RAGE.nii.gz mag 1 1
#fslroi sub-${SUBJECT}_ses-${SESSION}_MP2RAGE.nii.gz phase 0 1

# get MP2RAGE flat image
cd $PROGDIR/mprageconvert
matlab -nosplash -nodisplay -r "addpath $PROGDIR/mprageconvert/Nifti_tools; create_mp2rage_command('$HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/', 're.nii.gz ', 'im.nii.gz', 'sub-${SUBJECT}_ses-${SESSION}_MP2RAGE_DelRec1.json'); exit"

cd $HOMEDIR/data_BIDS/sub-${SUBJECT}/ses-${SESSION}/anat/
fslmaths MP2RAGE.nii.gz -add 0.5 MP2RAGEpos.nii.gz 
fslmaths mag -mul MP2RAGEpos magRAGE
 
fi

###############################################################################
# Before continuing, run structural analyses to get average T1w images
###############################################################################

###############################################################################
# Process with fmriprep
###############################################################################
if [ $PHASE == 0 ] || [ $PHASE == 3 ]; then

# create singularity container
if [ ! -e $HOMEDIR/fmriprep.simg ]; then
    SINGULARITY_TMPDIR=~/Data/tmp SINGULARITY_CACHEDIR=~/Data/tmp singularity build ~/Data/fmriprepversions/fmriprep20.1.1.simg docker://poldracklab/fmriprep:latest
    ln -s ~/Data/fmriprepversions/fmriprep20.1.1.simg $HOMEDIR/fmriprep.simg 
fi

rm -r $WORK
mkdir $WORK
cp /usr/local/freesurfer/license.txt $WD/

rm -r $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}

PYTHONPATH="" singularity run --bind /data:/data \
   --cleanenv $HOMEDIR/fmriprep.simg \
   $HOMEDIR/data_BIDS \
   $WD \
   participant \
   --ignore slicetiming \
   --longitudinal \
   --force-syn \
   --fs-license-file \
   $WD/license.txt \
   --nthreads $FMRIPREPPROCS \
   -w $WORK \
   --use-aroma \
   --write-graph \
   --skip_bids_validation \
   --output-spaces T1w MNI152NLin2009cAsym fsaverage6 \
   --participant-label ${SUBJECT} 

# clean up 
rm -r $WORK

fi # PHASE 3


###############################################################################
# Prepare surfaces, labels and masks
###############################################################################

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

# Import labels to subject 
LABELSDIR=$HOMEDIR/fmriprep/analysis/sub-$SUBJECT/label
mkdir $LABELSDIR

if [ ! -e $HOMEDIR/fmriprep/analysis/surf ]; then 
    mkdir $HOMEDIR/fmriprep/analysis/surf
fi
if [ ! -e $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf ]; then 
    mkdir $HOMEDIR/fmriprep/analysis/sub-$SUBJECT/surf
fi

for hemi in lh rh; do
  mri_surf2surf --srcsubject fsaverage \
        --srcsurfval $HOMEDIR/labels/fsaverage/${hemi}.somatomotor-mask.func.gii  \
        --trgsubject sub-${SUBJECT} \
        --trgsurfval $LABELSDIR/${hemi}.somatomotor-mask.func.gii \
        --hemi $hemi \
        --cortex

  mri_cor2label --i $LABELSDIR/${hemi}.somatomotor-mask.func.gii \
   --id 1 --l $LABELSDIR/${hemi}.somatomotor-mask.label \
   --surf sub-$SUBJECT $hemi white 
done

if [ -e $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_ses-${SESSION}_desc-brain_mask.nii.gz ]; then
    mri_vol2vol --mov $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_ses-${SESSION}_desc-brain_mask.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data_mean.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask.nii.gz --nearest
else 
    mri_vol2vol --mov $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/anat/sub-${SUBJECT}_desc-brain_mask.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data_mean.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask.nii.gz --nearest
fi

mri_label2vol --label $LABELSDIR/lh.somatomotor-mask.label --regheader T1.mgz --subject sub-${SUBJECT} --fill-ribbon --hemi lh \
--o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/lh.ribbon.nii.gz
mri_label2vol --label $LABELSDIR/rh.somatomotor-mask.label --regheader T1.mgz --subject sub-${SUBJECT} --fill-ribbon --hemi rh \
--o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/rh.ribbon.nii.gz

mri_vol2vol --mov $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/lh.ribbon.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data_mean.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/lh.ribbon.nii.gz --nearest
mri_vol2vol --mov $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/rh.ribbon.nii.gz --targ $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/data_mean.nii.gz --regheader --o $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/rh.ribbon.nii.gz --nearest

fslmaths $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/lh.ribbon.nii.gz -add \
$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/rh.ribbon.nii.gz -bin -kernel boxv 7 -dilF \
-mul $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask -bin  \
$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/mask
rm $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/*ribbon.nii.gz

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

if [ ! -e $LABELSDIR/${hemi}.somatomotor-mask.ds.label ]; then 
# get labels if they do not exist
for hemi in rh lh; do
    for surface in pial white inflated; do
        mri_surf2surf --hemi $hemi --srcsubject sub-${SUBJECT} --sval-xyz $surface --trgsubject fsaverage6 --trgicoorder 6 --trgsurfval  $SUBDIR/surf/${hemi}.${surface}.ds.mgh --tval-xyz $SUBJECTS_DIR/sub-${SUBJECT}/mri/T1.mgz 
        mris_convert --to-scanner $SUBDIR/surf/${hemi}.${surface}.ds.mgh $SUBDIR/surf/${hemi}.${surface}.ds.gii
    done

ln -s $SUBJECTS_DIR/fsaverage6/surf/${hemi}.sphere.reg $SUBJECTS_DIR/sub-${SUBJECT}/surf/${hemi}.sphere6.reg
ln -s $SUBDIR/surf/${hemi}.white.ds.mgh $SUBJECTS_DIR/sub-${SUBJECT}/surf/${hemi}.white.ds.mgh

# adapt mask
mri_convert $LABELSDIR/${hemi}.somatomotor-mask.func.gii $LABELSDIR/${hemi}.somatomotor-mask.func.mgh

mri_surf2surf --srcsubject sub-${SUBJECT} \
   --srcsurfval $LABELSDIR/${hemi}.somatomotor-mask.func.mgh \
   --srcsurfreg sphere.reg \
   --trgsubject sub-${SUBJECT} \
   --trgsurfval $LABELSDIR/${hemi}.somatomotor-mask.ds.func.mgh \
   --trgsurfreg sphere6.reg \
   --hemi $hemi 
    
mri_cor2label --i $LABELSDIR/${hemi}.somatomotor-mask.ds.func.mgh \
   --id 1 --l $LABELSDIR/${hemi}.somatomotor-mask.ds.label \
   --surf sub-$SUBJECT $hemi white.ds.mgh
done
rm $SUBDIR/surf/*.ds.mgh $LABELSDIR/*.somatomotor-mask.func.gii $LABELSDIR/*.somatomotor-mask*func.mgh $SUBJECTS_DIR/sub-${SUBJECT}/surf/${hemi}.white.ds.mgh
fi

fi #PHASE 4

if [ "$OVERWRITE" -eq 1 ]; then
    for run in $RUNS; 
    do    
          RUN_DIR=$SUBDIR/run$run
          rm -r $RUN_DIR
    done
fi

###############################################################################
# Run fMRI analysis
###############################################################################
if [ $PHASE == 0 ] || [ $PHASE == 5 ]; then

cd $SUBDIR

for run in $RUNS; 
do
  RUN_DIR=$SUBDIR/run$run

  if [ -e "$RUN_DIR/volume" ] && \
  [ ! -e "$RUN_DIR/volume/cope1.nii.gz" ] ; then
     rm -r $RUN_DIR/volume $RUN_DIR/surf*
  fi 
  if [ ! -e "$RUN_DIR/volume" ]; then

      FUNCVOL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
      FUNCSURFL=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_space-fsaverage6_hemi-L.func.gii # corresponds to fsaverage6 space
      FUNCSURFR=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_space-fsaverage6_hemi-R.func.gii

      mkdir $RUN_DIR
    
      #Smoothing       
      for HEMI in lh rh; do

      if [ "$HEMI" == "lh" ]; then
    	FUNCSURF=$FUNCSURFL
      else
    	FUNCSURF=$FUNCSURFR
      fi

      mris_convert $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.sphere.reg $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.sphere.reg.gii

      $HCPDIR/wb_shortcuts -freesurfer-resample-prep $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.white  \
      $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.pial \
      $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.sphere.reg \
      $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.sphere.reg.gii \
      $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.midthickness.surf.gii \
      $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.midthickness.164k_fs_LR.surf.gii \
      $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.sphere.reg.surf.gii

      # temporal filtering
      $HCPDIR/wb_command -metric-convert -to-nifti $FUNCSURF $RUN_DIR/fake_nifti.nii.gz
      fslmaths $RUN_DIR/fake_nifti.nii.gz -Tstd -bin $RUN_DIR/fake_mask.nii.gz
      FACTOR=`fslstats $RUN_DIR/fake_nifti.nii.gz -k $RUN_DIR/fake_mask.nii.gz -p 50`
      fslmaths $RUN_DIR/fake_nifti.nii.gz -mul $FACTOR $RUN_DIR/fake_nifti.nii.gz
      fslmaths $RUN_DIR/fake_nifti.nii.gz -Tmean $RUN_DIR/fake_mean.nii.gz
      fslmaths $RUN_DIR/fake_nifti.nii.gz -bptf `echo "0.5 * $TemporalFilter / $TR" | bc -l` -1 -add $RUN_DIR/fake_mean.nii.gz $RUN_DIR/fake_nifti.nii.gz
      
      $HCPDIR/wb_command -metric-convert -from-nifti $RUN_DIR/fake_nifti.nii.gz $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.midthickness.surf.gii $RUN_DIR/${HEMI}.filtered.func.gii
      rm $RUN_DIR/fake_*

      # avoid negative or zero values outside the cortex
      $HCPDIR/wb_command -metric-math "x*(x>=0) + 0*(x<0)" $RUN_DIR/${HEMI}.dilated.func.gii \
      -var x $RUN_DIR/${HEMI}.filtered.func.gii

      $HCPDIR/wb_command -metric-dilate $RUN_DIR/${HEMI}.dilated.func.gii $SUBJECTS_DIR/fsaverage6/surf/${HEMI}.midthickness.surf.gii 100 $RUN_DIR/${HEMI}.dilated.func.gii -nearest 

      done
    
      cut -d$'\t' -f $CONFOUND_INDICES $HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_desc-confounds_regressors.tsv | tail -n +2 > $RUN_DIR/mc_run-0${run}.csv
    
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
      
      NVOLS=`fslnvols $FUNCVOL`
      sed -i "s%@NVOLS%$NVOLS%" $FSF_FILE
      sed -i "s%@fMRI%$RUN_DIR/data_smooth.nii.gz%" $FSF_FILE
      sed -i "s%@STRUCTURAL_BRAIN%$STRUCT%" $FSF_FILE
      sed -i "s%@ANALYSIS%$RUN_DIR/volume/analysis.feat%" $FSF_FILE  
      sed -i "s%@CONFOUNDS%$RUN_DIR/mc_run-0${run}.csv%" $FSF_FILE        
      #sed -i "s%@FWHM%$FWHM%" $FSF_FILE        
      
      # now run it
      feat_model fMRI
      rm -r $RUN_DIR/surfL $RUN_DIR/surfR $RUN_DIR/volume 

      # surface	
      film_gls --rn=$RUN_DIR/surfL --sa --ms=15 --epith=$FWHM_SURF --in=$RUN_DIR/lh.dilated.func.gii --in2=$SUBJECTS_DIR/fsaverage6/surf/lh.midthickness.surf.gii --pd=fMRI.mat --con=fMRI.con --mode=surface  
      film_gls --rn=$RUN_DIR/surfR --sa --ms=15 --epith=$FWHM_SURF --in=$RUN_DIR/rh.dilated.func.gii --in2=$SUBJECTS_DIR/fsaverage6/surf/rh.midthickness.surf.gii --pd=fMRI.mat --con=fMRI.con --mode=surface  

      # volume, smooth and run
      fslmaths $FUNCVOL -s $SIGMA $RUN_DIR/data_smooth.nii.gz
      feat fMRI.fsf
      # film_gls --rn=$RUN_DIR/volume --sa --ms=$FWHM --in=$FUNCVOL--thr=1000 --pd=fMRI.mat --con=fMRI.con --mode=volume 

      mv $RUN_DIR/volume/analysis.feat/stats/*nii.gz $RUN_DIR/volume
      
      rm -r $RUN_DIR/volume/analysis.feat $RUN_DIR/*.dilated.func.gii $RUN_DIR/*.filtered.func.gii $RUN_DIR/data_smooth.nii.gz
  fi # volume exists
done 

fi # PHASE 5

###############################################################################
# Obtain single-trial parameter estimates
###############################################################################
# assumes PHASE5 already run

if [ $PHASE == 0 ] || [ $PHASE == 6 ]; then

# LEAST SQUARES ALL
if [ ! -e "$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/effects_LSS.nii.gz" ]; then

for run in $RUNS; 
do
  FUNC=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_space-T1w_desc-preproc_bold.nii.gz

  RUN_DIR=$SUBDIR/run$run
  MODEL_DIR=$RUN_DIR/model
  rm -r $MODEL_DIR
  mkdir $MODEL_DIR

  FSF_FILE=$MODEL_DIR/fMRI.fsf

  cd $MODEL_DIR
  NEVS_ORIG=$(($NTRIALS + 2))
  NEVS_REAL=$(($NEVS_ORIG * 2))
  NCONS_ORIG=1
  NCONS_REAL=1

  ~/Software/LeftHand/process/setup_GLM.sh $NEVS_ORIG $NEVS_REAL

  cat $SINGLE1_FSF_FILE fsfchunk $SINGLE2_FSF_FILE > $FSF_FILE
  rm fsfchunk
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/TRIAL* $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/FIXATION $MODEL_DIR/
  cp $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/run$run/STRETCH $MODEL_DIR/

  NVOLS=`fslnvols $FUNC`
  NVOXELS=`echo \`fslinfo $FUNC | awk 'NR>=3 && NR <= 6 {print $2F}' \` | awk '{ for(j=i=1; i<=NF; i++) j*=$i; print j; j=0 }'`

  sed -i "s%@NEVS_ORIG%$NEVS_ORIG%" $FSF_FILE
  sed -i "s%@NEVS_REAL%$NEVS_REAL%" $FSF_FILE
  sed -i "s%@NCONS_ORIG%$NCONS_ORIG%" $FSF_FILE
  sed -i "s%@NCONS_REAL%$NCONS_REAL%" $FSF_FILE
  sed -i "s%@NVOLS%$NVOLS%" $FSF_FILE
  sed -i "s%@NVOXELS%$NVOXELS%" $FSF_FILE
  sed -i "s%@fMRI%$FUNC%" $FSF_FILE
  sed -i "s%@ANALYSIS%$MODEL_DIR/analysis%" $FSF_FILE 

  sed -i "s%@EV$(($NTRIALS + 1))%$MODEL_DIR/FIXATION%" $FSF_FILE        
  sed -i "s%@EV$(($NTRIALS + 2))%$MODEL_DIR/STRETCH%" $FSF_FILE        
  sed -i "s%@CONFOUNDS%$RUN_DIR/mc_run-0${run}.csv%" $FSF_FILE        

  for EV in `seq $(($NTRIALS))`; 
  do  
     sed -i "s%@EV$EV%$MODEL_DIR/TRIAL${EV}%" $FSF_FILE        
  done

  feat $FSF_FILE
  
  # repeat to get prewhitened data
  #film_gls --in=$RUN_DIR/model/analysis.feat/filtered_func_data --rn=$RUN_DIR/model/analysis.feat/stats \
# --pd=$RUN_DIR/model/analysis.feat/design.mat --thr=1000.0 --sa -ms=5 --con=$RUN_DIR/model/analysis.feat/design.con --outputPWdata

  fslmerge -t $RUN_DIR/effects_LSA_$run `ls -v $RUN_DIR/model/analysis.feat/stats/${FX_FILE}*.nii.gz` 
  python $PROGDIR/compute_effects.py $RUN_DIR/model/analysis.feat/design.mat \
  $RUN_DIR/effects_LSA_$run.nii.gz \
  $RUN_DIR/model/analysis.feat/stats/res4d.nii.gz \
  $SUBDIR/mask.nii.gz \
  $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/sequences.csv \
  $run \
  $RUN_DIR/effects_LSS_${run}.nii.gz \
  $RUN_DIR/derivatives_LSS_${run}.nii.gz #--multivar

  rm -r $MODEL_DIR
done

# LEAST SQUARES SINGLE
if [ ]; then
#if [ ! -e "$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/effects.nii.gz" ]; then

for run in $RUNS; 
do

  FUNC=$HOMEDIR/fmriprep/fmriprep/sub-${SUBJECT}/ses-${SESSION}/func/sub-${SUBJECT}_ses-${SESSION}_task-sequence_run-0${run}_space-T1w_desc-preproc_bold.nii.gz

  RUN_DIR=$SUBDIR/run$run

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

  fslmerge -t $RUN_DIR/effects_LSS_$run `ls -v $RUN_DIR/model?/analysis.feat/stats/${FX_FILE}1.nii.gz $RUN_DIR/model??/analysis.feat/stats/${FX_FILE}1.nii.gz` 
  fslmerge -t $RUN_DIR/derivatives_LSS_$run `ls -v $RUN_DIR/model?/analysis.feat/stats/${FX_FILE}2.nii.gz $RUN_DIR/model??/analysis.feat/stats/${FX_FILE}2.nii.gz` 
  #rm -r $RUN_DIR/model*

done

fi # LSS extra computation

# Merge parameters
cd $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/

  fslmerge -t effects `ls -v run*/effects_LSS_?.nii.gz` 
  fslmerge -t derivatives `ls -v run*/derivatives_LSS_?.nii.gz` 

  fslmaths effects -Tmean -s 3 effects_mean
  fslmaths derivatives -Tmean -s 3 derivatives_mean
  
  echo "Total volumes : " `fslnvols effects.nii.gz`
  rm -r run*/model*
  rm run*/effects* run*/derivatives*


fi # check if effects already exist
fi # PHASE 7

###############################################################################
# Prepare engine for searchlight analysis
###############################################################################

# Find available fMRI runs, based on the data that are valid
cd $HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/
RUNS=`ls run*/volume/cope1.nii.gz  | cut -d'/' -f1 | cut -d'n' -f2`
RUNS=`echo $RUNS`
echo "Found following runs: $RUNS"

if [ $PHASE == 0 ] || [ $PHASE == 7 ]; then

# Prepare engine
for hemi in rh lh; do
    if [ ! -e "$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/surf/${hemi}-${RADIUS}-qe.gzipped.hdf5" ] && \
    [ -e "$HOMEDIR/fmriprep/analysis/sub-${SUBJECT}/ses-${SESSION}/effects.nii.gz" ]; then
    python $PROGDIR/prepare_queryengine.py \
    $SUBDIR \
    $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/sequences.csv \
    $hemi $RADIUS "$RUNS"
    fi
done

fi #PHASE 7

###############################################################################
# Searchlight analysis 
###############################################################################
if [ $PHASE == 0 ] || [ $PHASE == 8 ]; then

LABELSDIR=$HOMEDIR/fmriprep/analysis/sub-$SUBJECT/label

# Run searchlight
for hemi in rh lh; do
if [ ! -e "$SUBDIR/surf/$TESTDIR/$hemi.sl_within_spread_correlation_${RADIUS}.fsaverage.func.mgh" ]; then

    python $PROGDIR/surface_searchlight.py \
    $SUBDIR \
    $HOMEDIR/responses/sub-${SUBJECT}/ses-${SESSION}/sequences.csv  \
    $LABELSDIR/${hemi}.somatomotor-mask.label\
    $hemi $RADIUS $NPROC $TESTDIR "$RUNS"

    # Resample to average
    for metric in $METRICS; do
      mri_convert $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.gii $SUBDIR/surf/$TESTDIR/$hemi.sl_${metric}_${RADIUS}.func.mgh
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
fi
done

rm $SUBDIR/surf/$TESTDIR/*.mgh $SUBDIR/surf/$TESTDIR/qe.mask.nii.gz #$SUBDIR/surf/*gzipped.hdf5 

fi #PHASE 8

###############################################################################
# Plot surface results
###############################################################################

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

# plot individual maps
if [ $PHASE == 0 ] || [ $PHASE == 9 ]; then

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

fi #PHASE 9

exit 1

# average what we have and find rois
if [ $PHASE == 0 ] || [ $PHASE == 10 ]; then

LABELS="S_front_inf S_postcentral G_precentral S_front_sup Pole_occipital S_intrapariet_and_P_trans G_temporal_middle G_and_S_cingul-Mid-Ant"

# done every session, although it would be enough to do it for the last one

if [ "${SESSION}" -eq 7 ] ; then
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

done # metric
fi

#rm $LABELSDIR/*label*
fi #PHASE 11

###############################################
# Extract roi data
###############################################
# average what we have
if [ $PHASE == 0 ] || [ $PHASE == 11 ]; then
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


