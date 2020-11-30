#!/bin/sh 
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

#WD=/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/
#IMAGES=`ls mwp1rsub* | grep -v 3201 | grep -v 3105` # longitudinal

#WD=/home/benjamin.garzon/Data/LeftHand/Lund1/cat12cross/
#IMAGES=`ls mwp1sub* | grep -v 3105.1` # cross-sectional

#WD=/home/benjamin.garzon/Data/LeftHand/Lund1/cat12crossbias12_15/
#FWHM=12
#SURF_FWHM=15

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/cat12crossbias8_10/
FWHM=8
SURF_FWHM=10

if [ ]; then

SIGMA=`echo $FWHM/2.3548 | bc -l`

cd $WD/data/mri
IMAGES=`ls mwp1sub*` # cross-sectional  | grep -v 3202.3
echo $IMAGES

echo $IMAGES  |tr ' ' "\n" | cut -d'-' -f2 | cut -d'_' -f1 | cut -d'.' -f1,2 --output-delimiter=' ' > $WD/image_list.txt
#less $WD/image_list.txt

rm $WD/data_s${FWHM}.nii.gz
fslmerge -t $WD/data.nii.gz $IMAGES
fslmaths $WD/data.nii.gz -nan -s $SIGMA $WD/data_s${FWHM}.nii.gz
fslmaths $WD/data_s${FWHM}.nii.gz -Tmin -thr 0.01 -mul /home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/somatomotor-mask.MNI.VBM.nii.gz -bin $WD/mask
fslmaths $WD/data_s${FWHM}.nii.gz -Tmin -thr 0.01 -bin $WD/mask_whole_brain

#flirt -in $WD/data_s8.nii.gz -ref $WD/data_s8.nii.gz -out $WD/data_s8_4mm.nii.gz -interp trilinear -applyisoxfm 4  &
#flirt -in  $WD/mask_whole_brain -ref  $WD/mask_whole_brain -out  $WD/mask_whole_brain_4mm.nii.gz -interp nearestneighbour -applyisoxfm 4  &

#flirt -in $WD/data_s8.nii.gz -ref $WD/data_s8.nii.gz -out $WD/data_s8_3mm.nii.gz -interp trilinear -applyisoxfm 3  &
#flirt -in  $WD/mask_whole_brain -ref  $WD/mask_whole_brain -out  $WD/mask_whole_brain_3mm.nii.gz -interp nearestneighbour -applyisoxfm 3  &

#flirt -in $WD/data_s8.nii.gz -ref $WD/data_s8.nii.gz -out $WD/data_s8_2mm.nii.gz -interp trilinear -applyisoxfm 2  &
#flirt -in  $WD/mask_whole_brain -ref  $WD/mask_whole_brain -out  $WD/mask_whole_brain_2mm.nii.gz -interp nearestneighbour -applyisoxfm 2  &
fi

cd $WD/data/surf

IMAGES=`cat ../../image_list.txt | sed 's/ /./g'`

for hemi in lh rh; do

#SPHERE_SURF=/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces_32k/${hemi}.sphere.freesurfer.gii
#RESAMPLED_METRIC=${hemi}.thickness.resampled_32k.sub-${im}_MP2RAGE.gii
#SMOOTH_SURF=/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces_32k/${hemi}.central.freesurfer.gii

SPHERE_SURF=/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces/${hemi}.sphere.freesurfer.gii
SMOOTH_SURF=/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces/${hemi}.central.freesurfer.gii

for im in $IMAGES; do
echo $im
	   RESAMPLED_METRIC=${hemi}.thickness.resampled.sub-${im}_MP2RAGE.gii

           mri_convert ${hemi}.thickness.sub-${im}_MP2RAGE ${hemi}.thickness.sub-${im}_MP2RAGE.gii
           
           $HCPDIR/wb_command -metric-resample \
           ${hemi}.thickness.sub-${im}_MP2RAGE.gii \
           ${hemi}.sphere.sub-${im}_MP2RAGE.gii \
           $SPHERE_SURF \
           BARYCENTRIC \
	   $RESAMPLED_METRIC

           $HCPDIR/wb_command -metric-smoothing \
           $SMOOTH_SURF \
	   $RESAMPLED_METRIC \
           ${SURF_FWHM} \
           s${SURF_FWHM}.$RESAMPLED_METRIC
rm $RESAMPLED_METRIC ${hemi}.thickness.sub-${im}_MP2RAGE.gii
done
done

for hemi in lh rh; do
  echo $hemi
  rm $WD/${hemi}.thickness.txt
  for i in s${SURF_FWHM}.${hemi}.thickness.resampled.sub*MP2RAGE.gii; do
	 $HCPDIR/wb_command -metric-convert -to-nifti $i ${i}.nii.gz
         echo $i | cut -d'_' -f1 | cut -d'-' -f2 | cut -d'.' -f1,2 --output-delimiter=' ' >> $WD/${hemi}.thickness.txt
  done
fslmerge -t $WD/${hemi}.thickness.${SURF_FWHM}.nii.gz s${SURF_FWHM}.${hemi}.thickness.resampled.sub*nii.gz
fslmaths $WD/${hemi}.thickness.${SURF_FWHM}.nii.gz -Tstd -bin $WD/${hemi}.cortex.mask.nii.gz
rm s${SURF_FWHM}.${hemi}.*nii.gz
done

#freeview -f /home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces_32k/rh.inflated.freesurfer.gii:Overlay=GROUP_x_TRAINING+_p.func.gii

