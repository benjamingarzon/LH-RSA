#!/usr/bin/sh

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer
SUBJECTS="lue1101 lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"

HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
NSESSIONS=7
SMOOTH=10

for SUBJECT in $SUBJECTS; do
    for SESSION in `seq $NSESSIONS`; do
       for HEMI in rh lh; do            
           mris_convert -c $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.gii
    
           mris_convert $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii
    
           mris_convert $SUBJECTS_DIR/fsaverage/surf/${HEMI}.sphere.reg \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii

           mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/${HEMI}.white \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white.fsaverage.gii
           
           $HCPDIR/wb_command -metric-resample \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii \
           BARYCENTRIC \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.func.gii

           $HCPDIR/wb_command -metric-smoothing \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white.fsaverage.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.func.gii \
           $SMOOTH \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.$SMOOTH.func.gii
                          
           $HCPDIR/wb_command -metric-convert -to-nifti \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.$SMOOTH.func.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.nii.gz 
           
           rm $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.*gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white.fsaverage.gii

       done # hemi
    done
    
    # gather volumes 
    fslmerge -t $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/rh.thickness.nii.gz $SUBJECTS_DIR/sub-${SUBJECT}.?.long.sub-${SUBJECT}.base/surf/rh.thickness.nii.gz
    fslmerge -t $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/lh.thickness.nii.gz $SUBJECTS_DIR/sub-${SUBJECT}.?.long.sub-${SUBJECT}.base/surf/lh.thickness.nii.gz
    fslmaths $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/rh.thickness.nii.gz -Tstd $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/rh.thickness.std.nii.gz
    fslmaths $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/lh.thickness.nii.gz -Tstd $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/lh.thickness.std.nii.gz

done

# gather all
ls $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/rh.thickness.nii.gz | cut -d'-' -f2 | cut -d '.' --output-delimiter ' ' -f1,2 > $SUBJECTS_DIR/rh.thickness.txt
ls $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/lh.thickness.nii.gz | cut -d'-' -f2 | cut -d '.' --output-delimiter ' ' -f1,2 > $SUBJECTS_DIR/lh.thickness.txt

fslmerge -t $SUBJECTS_DIR/rh.thickness.nii.gz $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/rh.thickness.nii.gz 
fslmerge -t $SUBJECTS_DIR/lh.thickness.nii.gz $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/lh.thickness.nii.gz 

mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/lh.white \
$SUBJECTS_DIR/lh.white.gii

mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/rh.white \
$SUBJECTS_DIR/rh.white.gii

fslmaths $SUBJECTS_DIR/lh.thickness.nii.gz -Tmean $SUBJECTS_DIR/lh.thickness_mean.nii.gz
fslmaths $SUBJECTS_DIR/rh.thickness.nii.gz -Tmean $SUBJECTS_DIR/rh.thickness_mean.nii.gz

$HCPDIR/wb_command -metric-convert -from-nifti \
$SUBJECTS_DIR/lh.thickness_mean.nii.gz \
$SUBJECTS_DIR/lh.white.gii \
$SUBJECTS_DIR/lh.thickness_mean.func.gii

$HCPDIR/wb_command -metric-convert -from-nifti \
$SUBJECTS_DIR/rh.thickness_mean.nii.gz \
$SUBJECTS_DIR/rh.white.gii \
$SUBJECTS_DIR/rh.thickness_mean.func.gii

rm $SUBJECTS_DIR/?h.thickness_mean.nii.gz 

# mean of std
for GROUP in 1 2; do 
    fslmerge -t $SUBJECTS_DIR/rh.thickness.${GROUP}.std.nii.gz $SUBJECTS_DIR/sub-lue1${GROUP}*.base/surf/rh.thickness.std.nii.gz
    fslmerge -t $SUBJECTS_DIR/lh.thickness.${GROUP}.std.nii.gz $SUBJECTS_DIR/sub-lue1${GROUP}*.base/surf/lh.thickness.std.nii.gz
    fslmaths $SUBJECTS_DIR/rh.thickness.${GROUP}.std.nii.gz -Tmean $SUBJECTS_DIR/rh.thickness.${GROUP}.std.nii.gz
    fslmaths $SUBJECTS_DIR/lh.thickness.${GROUP}.std.nii.gz -Tmean $SUBJECTS_DIR/lh.thickness.${GROUP}.std.nii.gz
    
    $HCPDIR/wb_command -metric-convert -from-nifti \
    $SUBJECTS_DIR/lh.thickness.${GROUP}.std.nii.gz \
    $SUBJECTS_DIR/lh.white.gii \
    $SUBJECTS_DIR/lh.thickness.${GROUP}.std.func.gii
    
    $HCPDIR/wb_command -metric-convert -from-nifti \
    $SUBJECTS_DIR/rh.thickness.${GROUP}.std.nii.gz \
    $SUBJECTS_DIR/rh.white.gii \
    $SUBJECTS_DIR/rh.thickness.${GROUP}.std.func.gii

done 
rm $SUBJECTS_DIR/?h.thickness.*std.nii.gz $SUBJECTS_DIR/?h.thickness_mean.nii.gz #$SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/?h.thickness.nii.gz
rm $SUBJECTS_DIR/sub-*.base/surf/?h.thickness.?.std.nii.gz $SUBJECTS_DIR/?h.white.gii

cd $SUBJECTS_DIR/

mris_calc -o rh.thickness.std.diff.func.gii rh.thickness.1.std.func.gii sub rh.thickness.2.std.func.gii
mris_calc -o lh.thickness.std.diff.func.gii lh.thickness.1.std.func.gii sub lh.thickness.2.std.func.gii

freeview -f fsaverage/surf/lh.inflated:overlay=lh.thickness.std.diff.func.gii fsaverage/surf/rh.inflated:overlay=rh.thickness.std.diff.func.gii
#freeview -f fsaverage/surf/lh.inflated:overlay=lh.thickness_mean.func.gii fsaverage/surf/rh.inflated:overlay=rh.thickness_mean.func.gii

#freeview -f fsaverage/surf/lh.inflated:overlay=lh.thickness.1.std.func.gii fsaverage/surf/rh.inflated:overlay=rh.thickness.1.std.func.gii
#freeview -f fsaverage/surf/lh.inflated:overlay=lh.thickness.2.std.func.gii fsaverage/surf/rh.inflated:overlay=rh.thickness.2.std.func.gii

fslmaths $SUBJECTS_DIR/lh.thickness.nii.gz -Tmean -bin $SUBJECTS_DIR/lh.thickness.mask.nii.gz
fslmaths $SUBJECTS_DIR/rh.thickness.nii.gz -Tmean -bin $SUBJECTS_DIR/rh.thickness.mask.nii.gz

fslmaths $SUBJECTS_DIR/lh.thickness.nii.gz -Tmean -thr 3 -bin $SUBJECTS_DIR/lh.thickness.small.mask.nii.gz
fslmaths $SUBJECTS_DIR/rh.thickness.nii.gz -Tmean -thr 3 -bin $SUBJECTS_DIR/rh.thickness.smallmask.nii.gz

mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/lh.pial $SUBJECTS_DIR/lh.fsaverage.white.gii
mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/rh.pial $SUBJECTS_DIR/rh.fsaverage.white.gii

