#!/bin/sh
WB=/home/share/Software/HCP/workbench/bin_rh_linux64/wb_command
TEMPLATES_DIR=/home/share/Software/HCP/global/templates2/standard_mesh_atlases/resample_fsaverage

SOURCE_R=$TEMPLATES_DIR/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii
TARGET_R=$TEMPLATES_DIR/fsaverage_std_sphere.R.164k_fsavg_R.surf.gii
SOURCE_AREA_R=$TEMPLATES_DIR/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii
TARGET_AREA_R=$TEMPLATES_DIR/fsaverage.R.midthickness_va_avg.164k_fsavg_R.shape.gii

SOURCE_L=$TEMPLATES_DIR/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii
TARGET_L=$TEMPLATES_DIR/fsaverage_std_sphere.L.164k_fsavg_L.surf.gii
SOURCE_AREA_L=$TEMPLATES_DIR/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii
TARGET_AREA_L=$TEMPLATES_DIR/fsaverage.L.midthickness_va_avg.164k_fsavg_L.shape.gii

#SOURCE_L= /home/share/Software/HCP/global/templates/standard_mesh_atlases/fs_L/fsaverage.L.sphere.164k_fs_L.surf.gii #fs_LR.32k.L.white.surf.gii
#TARGET_L=/usr/local/freesurfer/subjects/fsaverage/surf/lh.white

TARGET_FILE=Icosahedron-162
#TARGET_FILE=Icosahedron-42

$WB -label-resample \
      ${TARGET_FILE}.32k.R.label.gii \
      $SOURCE_R \
      $TARGET_R \
      ADAP_BARY_AREA \
      ${TARGET_FILE}.fsaverage.R.label.gii \
      -area-metrics $SOURCE_AREA_R $TARGET_AREA_R 
     
$WB -label-resample \
      ${TARGET_FILE}.32k.L.label.gii \
      $SOURCE_L \
      $TARGET_L \
      ADAP_BARY_AREA \
      ${TARGET_FILE}.fsaverage.L.label.gii \
      -area-metrics $SOURCE_AREA_L $TARGET_AREA_L 

# convert to annotfile
mris_convert --annot ${TARGET_FILE}.fsaverage.R.label.gii /usr/local/freesurfer/subjects/fsaverage/surf/rh.white ./${TARGET_FILE}.fsaverage.R.annot
mris_convert --annot ${TARGET_FILE}.fsaverage.L.label.gii /usr/local/freesurfer/subjects/fsaverage/surf/lh.white ./${TARGET_FILE}.fsaverage.L.annot

# find the parcels that fall within the mask
SUBJECTS_DIR=/usr/local/freesurfer/subjects/
MASK_R=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/rh.somatomotor-mask.func.gii
MASK_L=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/lh.somatomotor-mask.func.gii
mri_segstats --annot fsaverage rh ./${TARGET_FILE}.fsaverage.R.label.gii --i $MASK_R --sum rh.somatomotor.labels
mri_segstats --annot fsaverage lh ./${TARGET_FILE}.fsaverage.L.label.gii --i $MASK_L --sum lh.somatomotor.labels
echo `awk '{if ($8 == 1.0000) print $5}' rh.somatomotor.labels` > rh.tessellation162.txt
echo `awk '{if ($8 == 1.0000) print $5}' lh.somatomotor.labels` > lh.tessellation162.txt

cat rh.tessellation162.txt | sed 's/label/rh.label/g' > MNI.tessellation162.txt
cat lh.tessellation162.txt | sed 's/label/lh.label/g' >> MNI.tessellation162.txt

cp *.tessellation.txt ~/Software/LeftHand/masks

TARGET_FILE=Icosahedron-162
freeview -f /usr/local/freesurfer/subjects/fsaverage/surf/rh.inflated:annot=${TARGET_FILE}.fsaverage.R.annot \
	    /usr/local/freesurfer/subjects/fsaverage/surf/lh.inflated:annot=${TARGET_FILE}.fsaverage.L.annot
	    
# MNI space


SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject
TARGET_FILE=Icosahedron-162
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/
LABELDIR=$SUBJECTS_DIR/fsaverage/label
PARCDIR=/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations/fs_LR_32
OUTPUTDIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/$TARGET_FILE
OUTPUTDIR2=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/cvs_avg35_inMNI152/$TARGET_FILE
mkdir $OUTPUTDIR
mkdir $OUTPUTDIR2

ln -s $PARCDIR/${TARGET_FILE}.fsaverage.R.annot $LABELDIR/rh.${TARGET_FILE}.annot
ln -s $PARCDIR/${TARGET_FILE}.fsaverage.L.annot $LABELDIR/lh.${TARGET_FILE}.annot

mri_annotation2label --subject fsaverage  --hemi rh --labelbase $OUTPUTDIR/rh.label --annotation ${TARGET_FILE} --ctab $OUTPUTDIR/rh.LUT.txt
mri_annotation2label --subject fsaverage  --hemi lh --labelbase $OUTPUTDIR/lh.label --annotation ${TARGET_FILE} --ctab $OUTPUTDIR/lh.LUT.txt

for hemi in rh lh; do
for i in  $(seq -f "%03g" 1 152); do
echo $i
if [ -f "$OUTPUTDIR2/${hemi}.label-${i}.label" ]; then 
  echo Label exists
  continue
fi
mri_label2label --srclabel $OUTPUTDIR/${hemi}.label-${i}.label \
--srcsubject fsaverage \
--trgsubject cvs_avg35_inMNI152 \
--trglabel $OUTPUTDIR2/${hemi}.label-${i}.label \
--hemi ${hemi} --regmethod surface

#mri_label2vol --label $OUTPUTDIR2/${hemi}.label-${i}.label \
#--regheader T1.mgz --subject cvs_avg35_inMNI152 --fill-ribbon --hemi $hemi \
#--o $OUTPUTDIR2/${hemi}.label-${i}.nii.gz

#mri_cor2label --i $OUTPUTDIR2/${hemi}.label-${i}.nii.gz \
#--id 1 --l $OUTPUTDIR2/${hemi}.label-${i}.vol.label 
done

mris_label2annot --s cvs_avg35_inMNI152 \
   --h $hemi \
   --ctab $OUTPUTDIR/LUT.padded.txt \
   --a $TARGET_FILE \
   --ldir $OUTPUTDIR2 \
   --sd $SUBJECTS_DIR 

done

mri_aparc2aseg --s cvs_avg35_inMNI152 --o rh.${TARGET_FILE}.aparc.nii.gz --annot $TARGET_FILE --rh --new-ribbon --rip-unknown --annot-table $OUTPUTDIR/LUT.padded.txt
mri_aparc2aseg --s cvs_avg35_inMNI152 --o lh.${TARGET_FILE}.aparc.nii.gz --annot $TARGET_FILE --lh --new-ribbon --rip-unknown --annot-table $OUTPUTDIR/LUT.padded.txt

#dilate slighly the labels (modal dilation)
TARGET_FILE=Icosahedron-162

mri_extract_label $SUBJECTS_DIR/cvs_avg35_inMNI152/mri/ribbon.mgz 42 rh.ribbon.nii.gz
mri_extract_label $SUBJECTS_DIR/cvs_avg35_inMNI152/mri/ribbon.mgz 3 lh.ribbon.nii.gz
fslmaths rh.${TARGET_FILE}.aparc.nii.gz -sub 2000 -thr 0 -mas rh.ribbon.nii.gz rh.${TARGET_FILE}.aparc.clean.nii.gz
fslmaths lh.${TARGET_FILE}.aparc.nii.gz -sub 1000 -thr 0 -mas lh.ribbon.nii.gz lh.${TARGET_FILE}.aparc.clean.nii.gz
fslmaths rh.${TARGET_FILE}.aparc.clean.nii.gz -dilD -kernel 3D rh.${TARGET_FILE}.aparc.dil.nii.gz
fslmaths lh.${TARGET_FILE}.aparc.clean.nii.gz -dilD -kernel 3D lh.${TARGET_FILE}.aparc.dil.nii.gz

fslmaths lh.${TARGET_FILE}.aparc.dil.nii.gz -add 152 -mas lh.${TARGET_FILE}.aparc.dil.nii.gz -add rh.${TARGET_FILE}.aparc.dil.nii.gz ${TARGET_FILE}.aparc.nii.gz 
mri_vol2vol --mov ${TARGET_FILE}.aparc.nii.gz --targ /home/benjamin.garzon/Data/LeftHand/Lund1/cat12crossbias8_10/data.nii.gz \
--o ${TARGET_FILE}.aparc.cat.nii.gz  --regheader --nearest

cat $OUTPUTDIR/LUT.padded.txt | awk -F ' ' '{print $1, $2, $3, $4, $5, $6}' | sed 's/label-/rh.label-/g' > $OUTPUTDIR2/LUT.txt
cat $OUTPUTDIR/LUT.padded.txt | awk -F ' ' '{print $1+152, $2, $3, $4, $5, $6}' | sed 's/label-/lh.label-/g' |sed 1d >> $OUTPUTDIR2/LUT.txt


