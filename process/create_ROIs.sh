#!/bin/sh
clear 
HOME=/home/xgarzb@GU.GU.SE/
WD=/data/lv0/MotorSkill
HCPDIR=/data/lv0/Software/workbench/bin_rh_linux64/
SUBJECTS_DIR=$WD/labels/subject
FS_DIR=/usr/local/freesurfer/7.1.1-1
MYSUBJECT=fsaverage
LABELSDIR=$WD/labels/$MYSUBJECT
PARC_DIR=/home/xgarzb@GU.GU.SE/Data/LeftHand/Lund1/labels/GlasserParc
RADIUS=25

SURF_L_MAP=$WD/fmriprep/analysis/higherlevel/Trained_Untrained/surfL/tests/quadratic/INTERCEPT_p.func.gii
SURF_R_MAP=$WD/fmriprep/analysis/higherlevel/Trained_Untrained/surfR/tests/quadratic/INTERCEPT_p.func.gii
SURF_L_AVERAGE=`echo ${SURF_L_MAP} | cut -f 1 -d '.'`.fsaverage.func.gii
SURF_R_AVERAGE=`echo ${SURF_R_MAP} | cut -f 1 -d '.'`.fsaverage.func.gii

if [ ]; then
# resample maps 
mri_surf2surf --srcsubject $MYSUBJECT --trgsubject fsaverage --sd $SUBJECTS_DIR \
--sval $SURF_R_MAP --tval $SURF_R_AVERAGE --hemi rh
mri_surf2surf --srcsubject $MYSUBJECT --trgsubject fsaverage --sd $SUBJECTS_DIR \
--sval $SURF_L_MAP --tval $SURF_L_AVERAGE --hemi lh

#freeview -f /data/lv0/MotorSkill/labels/subject/fsaverage/surf/rh.inflated:overlay=$SURF_R_AVERAGE:annot=$PARC_DIR/rh.HCP-MMP1.annot /data/lv0/MotorSkill/labels/subject/fsaverage/surf/lh.inflated:overlay=$SURF_L_AVERAGE:annot=$PARC_DIR/lh.HCP-MMP1.annot 
#exit 1

fi
cd $LABELSDIR
#primary sensorimotor
mri_mergelabels -i rh.R_1_ROI.label -i rh.R_3a_ROI.label -i rh.R_3b_ROI.label -i rh.R_4_ROI.label -o rh.wholePS.label
mri_mergelabels -i lh.L_1_ROI.label -i lh.L_3a_ROI.label -i lh.L_3b_ROI.label -i lh.L_4_ROI.label -o lh.wholePS.label

# create point
echo "#!ascii label  , from subject $MYSUBJECT vox2ras=TkReg" 1 > rh.point.label
echo "1" >> rh.point.label
echo "10669  34.83  -26.12  45.77 0.0" >> rh.point.label
echo "#!ascii label  , from subject $MYSUBJECT vox2ras=TkReg" 1 > lh.point.label
echo "1" >> lh.point.label
echo "15907  -34.95  -27.74  47.13 0.0" >> lh.point.label

for hemi in rh lh; do
mris_convert --label ${hemi}.point.label point $SUBJECTS_DIR/$MYSUBJECT/surf/${hemi}.pial ${hemi}.point.label.gii
mris_convert --label ${hemi}.wholePS.label wholePS $SUBJECTS_DIR/$MYSUBJECT/surf/${hemi}.pial ${hemi}.wholePS.label.gii
$HCPDIR/wb_command -gifti-label-to-roi ${hemi}.point.label.gii -name point ${hemi}.point.func.gii
$HCPDIR/wb_command -gifti-label-to-roi ${hemi}.wholePS.label.gii -name wholePS ${hemi}.wholePS.func.gii

mris_convert $SUBJECTS_DIR/$MYSUBJECT/surf/${hemi}.pial $SUBJECTS_DIR/$MYSUBJECT/surf/${hemi}.pial.shape.gii
rm ${hemi}.PS.func.gii ${hemi}.PS.label
#intersect with circular region
$HCPDIR/wb_command -metric-rois-from-extrema $SUBJECTS_DIR/$MYSUBJECT/surf/${hemi}.pial.shape.gii \
${hemi}.point.func.gii \
$RADIUS \
${hemi}.PS.func.gii \
-roi ${hemi}.wholePS.func.gii

mri_cor2label --i ${hemi}.PS.func.gii \
--id 1 --l $LABELSDIR/${hemi}.PS.label \
--surf $MYSUBJECT ${hemi} white --sd $SUBJECTS_DIR
done

freeview -f /data/lv0/MotorSkill/labels/subject/$MYSUBJECT/surf/rh.inflated:label=$LABELSDIR/rh.PS.label /data/lv0/MotorSkill/labels/subject/$MYSUBJECT/surf/lh.inflated:label=$LABELSDIR/lh.PS.label 
#exit 1
# superior parietal
mri_mergelabels -i rh.R_AIP_ROI.label -i rh.R_IP2_ROI.label -o rh.SPL.label
mri_mergelabels -i lh.L_AIP_ROI.label -i lh.L_IP2_ROI.label -o lh.SPL.label

# supplementary motor area
mri_mergelabels -i rh.R_6mp_ROI.label -i rh.R_6ma_ROI.label -i rh.R_SCEF_ROI.label -o rh.SMA.label 
mri_mergelabels -i lh.L_6mp_ROI.label -i lh.L_6ma_ROI.label -i lh.L_SCEF_ROI.label -o lh.SMA.label 

# premotor
rm rh.PM.label lh.PM.label rh.C1.label lh.C1.label rh.C2.label lh.C2.label
mri_mergelabels -i rh.R_FEF_ROI.label -i rh.R_6a_ROI.label -o rh.PM.label #cp rh.R_6a_ROI.label rh.PM.label
mri_mergelabels -i lh.L_FEF_ROI.label -i lh.L_6a_ROI.label -o lh.PM.label #cp lh.L_6a_ROI.label lh.PM.label

# a couple of control regions...
cp rh.R_PGi_ROI.label rh.C1.label 
cp lh.L_PGi_ROI.label lh.C1.label  

# create annot
rm /data/lv0/MotorSkill/labels/subject/$MYSUBJECT/label/rh.motor.rois.annot
rm /data/lv0/MotorSkill/labels/subject/$MYSUBJECT/label/lh.motor.rois.annot
#--l rh.PM.label \
mris_label2annot --s $MYSUBJECT --h rh --l rh.PS.label --l rh.SMA.label \
--l rh.SPL.label --l rh.PM.label --ctab colortable.txt --a motor.rois --sd /data/lv0/MotorSkill/labels/subject 
mris_label2annot --s $MYSUBJECT --h lh --l lh.PS.label --l lh.SMA.label \
--l lh.SPL.label --l lh.PM.label --ctab colortable.txt --a motor.rois --sd /data/lv0/MotorSkill/labels/subject 

for l in PS PM SMA SPL C1; do
  ln -s rh.${l}.label rh.R_${l}.label
  ln -s lh.${l}.label lh.L_${l}.label
done

# create mask with all the regions
mri_mergelabels -i rh.SMA.label -i rh.SPL.label -i rh.PM.label -i rh.PS.label -o rh.motor.rois.label
mri_mergelabels -i lh.SMA.label -i lh.SPL.label -i lh.PM.label -i lh.PS.label -o lh.motor.rois.label

for hemi in rh lh; do
mris_convert --label ${hemi}.motor.rois.label motor.rois $SUBJECTS_DIR/$MYSUBJECT/surf/${hemi}.pial ${hemi}.motor.rois.label.gii

$HCPDIR/wb_command -gifti-label-to-roi ${hemi}.motor.rois.label.gii -name motor.rois ${hemi}.motor.rois.func.gii
$HCPDIR/wb_command -metric-convert -to-nifti ${hemi}.motor.rois.func.gii ${hemi}.motor.rois.nii.gz

done
freeview -f /data/lv0/MotorSkill/labels/subject/$MYSUBJECT/surf/rh.inflated:annot=/data/lv0/MotorSkill/labels/subject/$MYSUBJECT/label/rh.motor.rois.annot /data/lv0/MotorSkill/labels/subject/$MYSUBJECT/surf/lh.inflated:annot=/data/lv0/MotorSkill/labels/subject/$MYSUBJECT/label/lh.motor.rois.annot
#less /usr/local/freesurfer/7.1.1-1/FreeSurferColorLUT.txt

