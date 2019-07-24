#!/bin/sh
# run freesurfers

EXPERT_FILE=/home/benjamin.garzon/Software/LeftHand/process/expert.opts

export SUBJECTS_DIR=$1
WD=$2
SUBJECT=$3

do_fs_cross(){
# check not done yet

SUB=$1
SESSION=$2
T1w=$3
T2w=$4
PD=$5
OUTDIR=$6
DIR=$7
#MASK=$DIR/template_brain.nii.gz
MASK=$SUBJECTS_DIR/${SUB}.template/mri/brainmask.mgz

if [ ! -f "$SUBJECTS_DIR//${SUB}.$SESSION/stats/lh.aparc.stats" ]; then 
    rm -r $SUBJECTS_DIR/${SUB}.$SESSION/
    echo "$SUBJECTS_DIR//${SUB}.$SESSION/ not done"
 
    # project the masks computed before       
    mri_vol2vol --mov $MASK \
    --targ $PD --lta-inv $OUTDIR/template.lta --nearest \
    --o $OUTDIR/brainmasknative.nii.gz

    mri_vol2vol --mov $T1w \
    --targ $MASK --lta $OUTDIR/template.lta --nearest \
    --o $DIR/ses-${SESSION}/anat/T1w_template.nii.gz
    
#    fslmaths $OUTDIR/brainmasknative.nii.gz -ero $OUTDIR/brainmasknative.nii.gz 
# clean the boundaries
    fslmaths $OUTDIR/brainmasknative.nii.gz -kernel boxv 7 -ero -bin $OUTDIR/min
    fslmaths $T1w -thr 0.7 -bin $OUTDIR/th
    fslmaths $OUTDIR/brainmasknative.nii.gz -bin -sub $OUTDIR/min -bin $OUTDIR/rim
    fslmaths $OUTDIR/rim -mul $OUTDIR/th -bin -sub 1 -mul -1 -mul $OUTDIR/brainmasknative.nii.gz $OUTDIR/brainmasknative.nii.gz  
    fslmaths $OUTDIR/brainmasknative.nii.gz -bin -kernel boxv 5 -fmean -thr 0.7 -bin $OUTDIR/brainmasknative.nii.gz
    fslmaths $T1w -mas  $OUTDIR/brainmasknative.nii.gz  $OUTDIR/T1wmasked.nii.gz 
    rm $OUTDIR/rim.nii.gz $OUTDIR/min.nii.gz $OUTDIR/th.nii.gz

    recon-all -autorecon1 -noskullstrip -s ${SUB}.$SESSION -i $OUTDIR/T1wmasked.nii.gz -hires \
    -expert $EXPERT_FILE

    cd $SUBJECTS_DIR/${SUB}.$SESSION/mri
    cp T1.mgz brainmask.auto.mgz
    ln -s brainmask.auto.mgz brainmask.mgz
    recon-all -autorecon2 -autorecon3 -s ${SUB}.$SESSION -T2 $T2w -T2pial -hires #-expert $EXPERT_FILE
    
else
    echo "$SUBJECTS_DIR/${SUB}.$SESSION/ cross already done"
fi

}

do_fs_long_base(){
SUB=$1
SESSIONS=$2
T2w=$3

# create base
if [ ! -f "$SUBJECTS_DIR/${SUB}.base/surf/lh.pial" ]; then 

    rm -r $SUBJECTS_DIR/${SUB}.base 
    
    recon-all -base ${SUB}.base -tp `echo $SESSIONS | sed 's/ / -tp /g'` -all -hires \
    -expert $EXPERT_FILE -T2 $T2w -T2pial 
else
   echo "$SUBJECTS_DIR/$SUB/ base already done"
fi

} 

do_fs_long(){
SUBDIR=$1
BASE=$2

# create longs
if [ ! -f "$SUBJECTS_DIR/${SUBDIR}.long.${BASE}/stats/lh.aparc.stats" ]; then 
  rm -r $SUBJECTS_DIR/${SUBDIR}.long.${BASE}
  recon-all -long ${SUBDIR} $BASE -all -T2pial -hires -expert $EXPERT_FILE
else
   echo "$SUBDIR long already done"
fi

}

##############################################################################

echo Doing subject $SUBJECT
ANATLIST=$WD/sub-$SUBJECT/ses-*/anat
NANAT=`echo $ANATLIST | wc -w`

STRUCTLIST=$WD/sub-$SUBJECT/ses-*/anat/magRAGE.nii.gz

# create template and do skull_stripping
NOBIAS=`echo $STRUCTLIST | sed 's/magRAGE.nii.gz/magRAGE_restore.nii.gz/g'`
TEMPLATES=`echo $STRUCTLIST | sed 's/magRAGE.nii.gz/mag_template.nii.gz/g'`
WEIGHTS=`echo $STRUCTLIST | sed 's/magRAGE.nii.gz/weights.nii.gz/g'`
LTAS=`echo $STRUCTLIST | sed 's/magRAGE.nii.gz/template.lta/g'`

if [ ! -f "$WD/sub-$SUBJECT/template_brain.nii.gz" ]; then
for STRUCT in $STRUCTLIST; do
    fast -B -v $STRUCT 
    rm $WD/sub-$SUBJECT/ses-*/anat/*pve*
    rm $WD/sub-$SUBJECT/ses-*/anat/*mixel*
    rm $WD/sub-$SUBJECT/ses-*/anat/*_seg.nii.gz
done

mri_robust_template --mov $NOBIAS --template $WD/sub-$SUBJECT/template.nii.gz --satit --lta $LTAS --mapmov $TEMPLATES --iscale --weights $WEIGHTS --maxit 30
bet $WD/sub-$SUBJECT/template.nii.gz $WD/sub-$SUBJECT/template_brain.nii.gz
rm -r $SUBJECTS_DIR/sub-${SUBJECT}.template
recon-all -autorecon1 -s sub-${SUBJECT}.template -i $WD/sub-$SUBJECT/template.nii.gz -hires 

fi

for DIR in $ANATLIST; do
    SESSION=`echo $DIR | cut -d'/' -f9 | cut -d '-' -f2` 
    
    T1w=$DIR/MP2RAGEpos.nii.gz 
    PD=$DIR/mag.nii.gz 
    T2w=$DIR/sub-${SUBJECT}_ses-${SESSION}_T2w.nii.gz 
    
    if [ -e "$T1w" ] && [ -e "$T2w" ] && [ -e "$PD" ] ; then
        # all available
        echo "Running cross-sectional"    
        do_fs_cross sub-${SUBJECT} $SESSION $T1w $T2w $PD $DIR $WD/sub-$SUBJECT/ &
    fi

done

# wait for all to finish
while [ `ls $SUBJECTS_DIR/sub-${SUBJECT}.?/stats/lh.aparc.stats | wc -w` -lt $NANAT ]; do
    echo "Waiting for all cross to finish"
    echo "sub-${SUBJECT}: Only `ls $SUBJECTS_DIR/sub-${SUBJECT}.?/stats/lh.aparc.stats | wc -w` out of $NANAT have finished"
    sleep 3600
done

if [ ! -e $WD/sub-$SUBJECT/template_T2w.nii.gz ]; then
  mri_robust_template --mov $WD/sub-$SUBJECT/ses-*/anat/sub-${SUBJECT}_ses-*_T2w.nii.gz --template $WD/sub-$SUBJECT/template_T2w.nii.gz --satit  --iscale --maxit 30
fi

# run longit
echo "Running base"    
do_fs_long_base sub-${SUBJECT} "$SUBJECTS_DIR/sub-${SUBJECT}.?" $WD/sub-$SUBJECT/template_T2w.nii.gz
ln -sf $SUBJECTS_DIR/sub-${SUBJECT}.base $SUBJECTS_DIR/sub-$SUBJECT 

fslmerge -t $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_all.nii.gz $WD/sub-$SUBJECT/ses-*/anat/T1w_template.nii.gz
fslmaths $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_all.nii.gz -log -Tmean -exp $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz
mri_vol2vol --mov $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz --targ $WD/sub-${SUBJECT}/template.nii.gz --nearest --regheader --out $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz  
fslmaths $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz -mul $WD/sub-${SUBJECT}/template.nii.gz $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w.nii.gz
rm $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_all.nii.gz $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz

# copy to anat dir
for DIR in $ANATLIST; do
    SESSION=`echo $DIR | cut -d'/' -f9 | cut -d '-' -f2` 
    ln -sf $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w.nii.gz $WD/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_ses-${SESSION}_T1w.nii.gz
done



cd $SUBJECTS_DIR
for SUBDIR in sub-$SUBJECT.?; do
        
        # all available
        echo "Running longitudinal"    
        do_fs_long ${SUBDIR} sub-${SUBJECT}.base &

done
