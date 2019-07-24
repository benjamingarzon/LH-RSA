#!/usr/bin/sh

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer
SUBJECTS="lue1101 lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
#SUBJECTS="lue1101" 
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
NSESSIONS=7


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
           
           $HCPDIR/wb_command -metric-resample \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii \
           BARYCENTRIC \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.func.gii
                          
           $HCPDIR/wb_command -metric-convert -to-nifti \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.func.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.nii.gz 
           
           rm $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness.*gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii

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

rm $SUBJECTS_DIR/rh.thickness.nii.gz $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/rh.thickness.nii.gz
rm $SUBJECTS_DIR/rh.thickness.nii.gz $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/rh.thickness.nii.gz

# mean of std
fslmerge -t $SUBJECTS_DIR/rh.thickness.std.nii.gz $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/rh.thickness.std.nii.gz
fslmerge -t $SUBJECTS_DIR/lh.thickness.std.nii.gz $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/lh.thickness.std.nii.gz
fslmaths $SUBJECTS_DIR/rh.thickness.std.nii.gz -Tmean $SUBJECTS_DIR/rh.thickness.std.nii.gz
fslmaths $SUBJECTS_DIR/lh.thickness.std.nii.gz -Tmean $SUBJECTS_DIR/lh.thickness.std.nii.gz
rm $SUBJECTS_DIR/sub-*.base/surf/lh.thickness.std.nii.gz
rm $SUBJECTS_DIR/sub-*.base/surf/rh.thickness.std.nii.gz

$HCPDIR/wb_command -metric-convert -from-nifti \
$SUBJECTS_DIR/lh.thickness.std.nii.gz \
/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/fsaverage/surf/lh.white \
$SUBJECTS_DIR/lh.thickness.std.gii

$HCPDIR/wb_command -metric-convert -from-nifti \
$SUBJECTS_DIR/rh.thickness.std.nii.gz \
/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/fsaverage/surf/rh.white \
$SUBJECTS_DIR/rh.thickness.std.gii

rm $SUBJECTS_DIR/?h.thickness.std.nii.gz