#!/bin/sh
# run freesurfers and prepare for VBM
##
#generated images
#sub-*/T1w_template.nii.gz: template of T1w images
#sub-*/T1w_template_brain.nii.gz

#sub-*/T2w_template.nii.gz: template of T2w images
#sub-*/T2w_template_brain.nii.gz

#sub-*/T2wtoT1w_template.nii.gz: T2w template registered to T1w template
#sub-*/T2wtoT1w_template_brain.nii.gz: Skull-stripped T2w template registered to T1w template
#sub-*/T2wtoT1w_template.lta: trasnformation of T2w template registered to T1w template

#sub-*/brainmask.nii.gz: T1w brain mask
#sub-*/ses-*/brainmasknative.nii.gz: T1w brain mask in original timepoint space

#sub-*/ses-*/T1wtotemplate.nii.gz: individual T1w in T1w space
#sub-*/ses-*/T1wtotemplate_brain.nii.gz: individual T1w in T1w space, skull-stripped
#sub-*/ses-*/T2wtoT1w.nii.gz: individual T2w in original timepoint space
#sub-*/ses-*/T2wtoT1w_template.nii.gz: individual T2w in T1w template space
#sub-*/ses-*/T2wtoT1w_template_brain.nii.gz: individual T2w in T1w template space, skull-stripped

# magRAGE mix of PD and MP2RAGE weight

export SUBJECTS_DIR=$1
WD=$2 # BIDS dir
STRUCT_DIR=$3 # BIDS dir
SUBJECT=$4
EXPERT_FILE=$5
OVERWRITE=$6
export DOFS=$7

REGCOST="NMI"

RECONALL=~/Software/LeftHand/process/recon-allT2

if [ "$OVERWRITE" -eq 1 ]; then
rm $STRUCT_DIR/sub-$SUBJECT/*template* $WD/sub-$SUBJECT/ses-*/anat/*restore*
rm -r $SUBJECTS_DIR/sub-$SUBJECT*
fi

do_fs_cross(){
# check not done yet

SUB=$1
SESSION=$2
T1w=$3
T2w=$4
PD=$5 # not used
DIR=$6

OUTDIR=$DIR/ses-$SESSION


MASK=$DIR/brainmask.nii.gz

if [ ! -f "$SUBJECTS_DIR/${SUB}.$SESSION/stats/lh.aparc.stats" ]; then 
    rm -r $SUBJECTS_DIR/${SUB}.$SESSION/
    echo "$SUBJECTS_DIR/${SUB}.$SESSION/ not done"
 
    # project the masks computed before       
    mri_vol2vol --mov $MASK \
    --targ $T1w --lta-inv $OUTDIR/T1wtotemplate.lta --nearest \
    --o $OUTDIR/brainmasknative.nii.gz    

    mri_vol2vol --mov $T1w \
    --targ $MASK --lta $OUTDIR/T1wtotemplate.lta --interp cubic \
    --o $OUTDIR/T1wtotemplate.nii.gz  # For vbm

    fslmaths $OUTDIR/T1wtotemplate.nii.gz -mas $MASK \
    $OUTDIR/T1wtotemplate_brain.nii.gz # for vbm
    
    mri_concatenate_lta $OUTDIR/T2wtotemplate.lta  $DIR/T2wtoT1w_template.lta $OUTDIR/compos1.lta     

    mri_vol2vol --mov $T2w \
    --targ $MASK --lta $OUTDIR/compos1.lta \
    --o $OUTDIR/T2wtoT1w_template.nii.gz --interp cubic  # for vbm
    
    fslmaths $OUTDIR/T2wtoT1w_template.nii.gz -mas $MASK \
    $OUTDIR/T2wtoT1w_template_brain.nii.gz # for vbm    
    
    mri_concatenate_lta -invert2 $OUTDIR/compos1.lta $OUTDIR/T1wtotemplate.lta $OUTDIR/compos2.lta

    mri_vol2vol --mov $T2w \
    --targ $T1w --lta $OUTDIR/compos2.lta \
    --o $OUTDIR/T2wtoT1w.nii.gz --no-resample    # for freesurfer
        
    # compute myelin maps
    fslmaths $OUTDIR/T1wtotemplate_brain.nii.gz -div $OUTDIR/T2wtoT1w_template_brain.nii.gz -mas $MASK $OUTDIR/myelin.nii.gz

    if [ "$DOFS" -eq 1 ]; then
    # clean the boundaries
    fslmaths $OUTDIR/brainmasknative.nii.gz -kernel boxv 7 -ero -bin $OUTDIR/min
    fslmaths $T1w -thr 0.7 -bin $OUTDIR/th
    fslmaths $OUTDIR/brainmasknative.nii.gz -bin -sub $OUTDIR/min -bin $OUTDIR/rim
    fslmaths $OUTDIR/rim -mul $OUTDIR/th -bin -sub 1 -mul -1 -mul $OUTDIR/brainmasknative.nii.gz $OUTDIR/brainmasknative.nii.gz  
    fslmaths $OUTDIR/brainmasknative.nii.gz -bin -kernel boxv 5 -fmean -thr 0.7 -bin $OUTDIR/brainmasknative.nii.gz
    rm $OUTDIR/rim.nii.gz $OUTDIR/min.nii.gz $OUTDIR/th.nii.gz 
    
    fslmaths $T1w -mas $OUTDIR/brainmasknative.nii.gz  $OUTDIR/T1wmasked.nii.gz 
    
    $RECONALL -autorecon1 -noskullstrip -s ${SUB}.$SESSION -i $OUTDIR/T1wmasked.nii.gz -hires \
    -expert $EXPERT_FILE

    cd $SUBJECTS_DIR/${SUB}.$SESSION/mri
    cp T1.mgz brainmask.auto.mgz
    ln -s brainmask.auto.mgz brainmask.mgz
    $RECONALL -autorecon2 -autorecon3 -s ${SUB}.$SESSION -T2 $OUTDIR/T2wtoT1w.nii.gz -T2pial -hires #-expert $EXPERT_FILE

    else 
       exit 1
    fi
    
    
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
    
    $RECONALL -base ${SUB}.base -tp `echo $SESSIONS | sed 's/ / -tp /g'` -all -hires \
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
  $RECONALL -long ${SUBDIR} $BASE -all -hires -expert $EXPERT_FILE #-T2pial
  echo "$SUBDIR is finished!"
  rm -r $SUBJECTS_DIR/${SUBDIR}
else
   echo "$SUBDIR long already done"
fi

}


##############################################################################

echo Doing subject $SUBJECT
ANATLIST=$WD/sub-$SUBJECT/ses-*/anat
NANAT=`echo $ANATLIST | wc -w`
SESLIST=`echo $ANATLIST | sed 's/anat//g' | sed "s%$WD%$STRUCT_DIR%g"`

T1wLIST=$WD/sub-$SUBJECT/ses-*/anat/magRAGE.nii.gz
T2wLIST=$WD/sub-$SUBJECT/ses-*/anat/sub-*_ses-*_T2w.nii.gz 

mkdir -p $STRUCT_DIR
mkdir -p $STRUCT_DIR/sub-$SUBJECT
# create structural directories
for SES in $SESLIST; do
  mkdir -p $SES
done

##################################################################
# create templates and do skull_stripping
##################################################################
T1wtoTEMPLATE=`echo $T1wLIST | sed 's%anat/magRAGE.nii.gz%T1wtotemplate.nii.gz%g'| sed "s%$WD%$STRUCT_DIR%g"`
T2wtoTEMPLATE=`echo $T1wLIST | sed 's%anat/magRAGE.nii.gz%T2wtotemplate.nii.gz%g'| sed "s%$WD%$STRUCT_DIR%g"`

T1wNOBIAS=`echo $T1wLIST | sed 's/.nii.gz/_restore.nii.gz/g'`
T2wNOBIAS=`echo $T2wLIST | sed 's/.nii.gz/_restore.nii.gz/g'`

WEIGHTS=`echo $T1wLIST | sed 's%anat/magRAGE.nii.gz%weights.nii.gz%g' | sed "s%$WD%$STRUCT_DIR%g"`
T1wLTAS=`echo $T1wLIST | sed 's%anat/magRAGE.nii.gz%T1wtotemplate.lta%g' | sed "s%$WD%$STRUCT_DIR%g"`
T2wLTAS=`echo $T1wLIST | sed 's%anat/magRAGE.nii.gz%T2wtotemplate.lta%g' | sed "s%$WD%$STRUCT_DIR%g"`


# build T1w template
if [ ! -f "$STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz" ]; then
for STRUCT in $T1wLIST; do
    fast -B -v $STRUCT &    
done

while [ `ls $WD/sub-$SUBJECT/ses-*/anat/magRAGE_restore.nii.gz | wc -l` -lt $NANAT ]; do
    echo "Waiting for T1w segmentations to finish"
    sleep 500
done

mri_robust_template --mov $T1wNOBIAS --template $STRUCT_DIR/sub-$SUBJECT/T1w_template.nii.gz --satit --lta $T1wLTAS --mapmov $T1wtoTEMPLATE --iscale --weights $WEIGHTS --maxit 30
rm -r $SUBJECTS_DIR/sub-${SUBJECT}.template
$RECONALL -autorecon1 -s sub-${SUBJECT}.template -i $STRUCT_DIR/sub-$SUBJECT/T1w_template.nii.gz -hires 
rm $WD/sub-$SUBJECT/ses-*/anat/*pve* $WD/sub-$SUBJECT/ses-*/anat/*mixel* $WD/sub-$SUBJECT/ses-*/anat/*_seg.nii.gz

mri_vol2vol --mov $SUBJECTS_DIR/sub-${SUBJECT}.template/mri/brainmask.mgz \
--targ $STRUCT_DIR/sub-$SUBJECT/T1w_template.nii.gz \
--out $STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz --regheader

fi

# build T2w template
if [ ! -e $STRUCT_DIR/sub-$SUBJECT/T2w_template.nii.gz ]; then

for STRUCT in $T2wLIST; do
    fast -B -v -t 2 $STRUCT &     
done

while [ `ls $WD/sub-$SUBJECT/ses-*/anat/*T2w_restore.nii.gz | wc -l` -lt $NANAT ]; do
    echo "Waiting for T2w segmentations to finish"
    sleep 500
done

rm $WD/sub-$SUBJECT/ses-*/anat/*pve* $WD/sub-$SUBJECT/ses-*/anat/*mixel* $WD/sub-$SUBJECT/ses-*/anat/*_seg.nii.gz
mri_robust_template --mov $T2wNOBIAS --template $STRUCT_DIR/sub-$SUBJECT/T2w_template.nii.gz --satit --lta $T2wLTAS  --mapmov $T2wtoTEMPLATE --iscale --maxit 30

bet $STRUCT_DIR/sub-$SUBJECT/T2w_template.nii.gz $STRUCT_DIR/sub-$SUBJECT/T2w_template_brain.nii.gz -R

# register T2w template to T1w template
mri_robust_register --mov $STRUCT_DIR/sub-$SUBJECT/T2w_template_brain.nii.gz \
--dst $STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz --satit --iscale \
--lta $STRUCT_DIR/sub-$SUBJECT/T2wtoT1w_template.lta --cost $REGCOST

mri_vol2vol --mov $STRUCT_DIR/sub-$SUBJECT/T2w_template_brain.nii.gz \
--targ $STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz \
--reg $STRUCT_DIR/sub-$SUBJECT/T2wtoT1w_template.lta \
--nearest --o $STRUCT_DIR/sub-$SUBJECT/T2wtoT1w_template_brain.nii.gz 

mri_vol2vol --mov $STRUCT_DIR/sub-$SUBJECT/T2w_template.nii.gz \
--targ $STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz \
--reg $STRUCT_DIR/sub-$SUBJECT/T2wtoT1w_template.lta \
--cubic --o $STRUCT_DIR/sub-$SUBJECT/T2wtoT1w_template.nii.gz 

# clean T1w further with skull-stripped T2w
fslmaths $STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz -mas $STRUCT_DIR/sub-$SUBJECT/T2wtoT1w_template_brain.nii.gz \
$STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz

fi

if [ ! -e $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz ]; then

# create a mask and mean volume
fslmaths $STRUCT_DIR/sub-$SUBJECT/T1w_template_brain.nii.gz -bin $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz
fslmerge -t $STRUCT_DIR/sub-$SUBJECT/T1wall.nii.gz $STRUCT_DIR/sub-$SUBJECT/ses-*/T1wtotemplate.nii.gz
fslmaths $STRUCT_DIR/sub-$SUBJECT/T1wall.nii.gz -Tmean $STRUCT_DIR/sub-$SUBJECT/T1wmean.nii.gz
rm $STRUCT_DIR/sub-$SUBJECT/T1wall.nii.gz

# clean the boundaries
fslmaths $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz -kernel boxv 7 -ero -bin $STRUCT_DIR/sub-$SUBJECT/min
fslmaths $STRUCT_DIR/sub-$SUBJECT/T1wmean.nii.gz -thr 0.7 -bin $STRUCT_DIR/sub-$SUBJECT/th
fslmaths $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz -bin -sub $STRUCT_DIR/sub-$SUBJECT/min -bin $STRUCT_DIR/sub-$SUBJECT/rim
fslmaths $STRUCT_DIR/sub-$SUBJECT/rim -mul $STRUCT_DIR/sub-$SUBJECT/th -bin \ 
-sub 1 -mul -1 -mul $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz  
fslmaths $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz -bin -kernel \
boxv 5 -fmean -thr 0.7 -bin $STRUCT_DIR/sub-$SUBJECT/brainmask.nii.gz
rm $STRUCT_DIR/sub-$SUBJECT/rim.nii.gz $STRUCT_DIR/sub-$SUBJECT/min.nii.gz $STRUCT_DIR/sub-$SUBJECT/th.nii.gz 
fi

##################################################################
# run cross-sectional reconstruction
##################################################################

for ANATDIR in $ANATLIST; do
    SESSION=`echo $ANATDIR | cut -d'/' -f9 | cut -d '-' -f2` 

    T1w=$ANATDIR/MP2RAGEpos.nii.gz 
    PD=$ANATDIR/mag.nii.gz 
    T2w=$ANATDIR/sub-${SUBJECT}_ses-${SESSION}_T2w.nii.gz 
    
    if [ -e "$T1w" ] && [ -e "$T2w" ] && [ -e "$PD" ] ; then
        # all available
        echo "Running cross-sectional"    
        do_fs_cross sub-${SUBJECT} $SESSION $T1w $T2w $PD $STRUCT_DIR/sub-$SUBJECT/ &
    fi

done

if [ ]; then
# wait for all to finish
while [ `ls $SUBJECTS_DIR/sub-${SUBJECT}.?/stats/lh.aparc.stats | wc -w` -lt $NANAT ]; do
    echo "Waiting for all cross to finish"
    echo "sub-${SUBJECT}: Only `ls $SUBJECTS_DIR/sub-${SUBJECT}.?/stats/lh.aparc.stats | wc -w` out of $NANAT have finished"
    sleep 500
done
fi

echo "Cross-sectionals done!"
##################################################################
# run base reconstruction
##################################################################
echo "Running base"    
do_fs_long_base sub-${SUBJECT} "$SUBJECTS_DIR/sub-${SUBJECT}.?" $STRUCT_DIR/sub-$SUBJECT/T2wtoT1w_template_brain.nii.gz
ln -sf $SUBJECTS_DIR/sub-${SUBJECT}.base $SUBJECTS_DIR/sub-$SUBJECT 

fslmerge -t $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_all.nii.gz $WD/sub-$SUBJECT/ses-*/anat/T1w_template.nii.gz
fslmaths $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_all.nii.gz -log -Tmean -exp $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz
mri_vol2vol --mov $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz --targ $WD/sub-${SUBJECT}/T1w_template.nii.gz --nearest --regheader --out $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz  
fslmaths $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz -mul $WD/sub-${SUBJECT}/T1w_template.nii.gz $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w.nii.gz
rm $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_all.nii.gz $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w_mean.nii.gz

# copy to anat dir
for DIR in $ANATLIST; do
    SESSION=`echo $DIR | cut -d'/' -f9 | cut -d '-' -f2` 
    ln -sf $WD/sub-${SUBJECT}/sub-${SUBJECT}_T1w.nii.gz $WD/sub-${SUBJECT}/ses-${SESSION}/anat/sub-${SUBJECT}_ses-${SESSION}_T1w.nii.gz
done

##################################################################
# run longitudinal reconstruction
##################################################################

cd $SUBJECTS_DIR
for SUBDIR in sub-$SUBJECT.?; do
        
        # all available
        echo "Running longitudinal"    
        do_fs_long ${SUBDIR} sub-${SUBJECT}.base &

done
