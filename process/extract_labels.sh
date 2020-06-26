#!/bin/sh

#labels of interest
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/
SUBJECTS_DIR=/usr/local/freesurfer/subjects/

LABELSDIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage
mkdir $LABELSDIR
cd $LABELSDIR

mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi lh --annotation aparc.a2009s --subject fsaverage
mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR  --hemi rh --annotation aparc.a2009s --subject fsaverage


LABELS="S_front_inf G_postcentral S_central S_postcentral G_precentral S_front_sup Pole_occipital S_intrapariet_and_P_trans G_temporal_middle G_and_S_cingul-Mid-Ant S_intrapariet_and_P_trans G_front_sup S_precentral-sup-part G_parietal_sup G_and_S_paracentral S_front_sup S_cingul-Marginalis G_and_S_cingul-Mid-Ant G_and_S_cingul-Mid-Post G_precuneus G_front_middle S_precentral-inf-part"

for hemi in lh rh; do
 
  for label in $LABELS; do
    echo $label
    mris_convert --label ${hemi}.${label}.label ${label}\
        $SUBJECTS_DIR/fsaverage/surf/${hemi}.white \
        ${hemi}.${label}.label.gii

    $HCPDIR/wb_command -gifti-label-to-roi ${hemi}.${label}.label.gii -name ${label} ${hemi}.${label}.func.gii
        
    $HCPDIR/wb_command -metric-convert -to-nifti \
    ${hemi}.${label}.func.gii \
    ${hemi}.${label}.nii.gz

  done
done

# gather a few in a mask
LABELS="S_precentral-sup-part G_postcentral S_central G_precentral S_postcentral S_intrapariet_and_P_trans G_front_sup G_parietal_sup G_and_S_paracentral S_front_sup S_cingul-Marginalis G_and_S_cingul-Mid-Ant G_and_S_cingul-Mid-Post G_precuneus G_front_middle S_precentral-inf-part"

rm ?h.mask.func.gii 
for hemi in lh rh; do
 
  for label in $LABELS; do
    echo $label

    if [ -e ${hemi}.mask.func.gii ]; then

      $HCPDIR/wb_command -metric-math "(x+y)" ${hemi}.mask.func.gii \
        -var x ${hemi}.mask.func.gii \
        -var y ${hemi}.${label}.func.gii  
    else 
       cp ${hemi}.${label}.func.gii ${hemi}.mask.func.gii
    fi
  done
    mv ${hemi}.mask.func.gii ${hemi}.somatomotor-mask.func.gii
    $HCPDIR/wb_command -metric-convert -to-nifti \
    ${hemi}.somatomotor-mask.func.gii \
    ${hemi}.mask.nii.gz

    # resample to fsaverage6 space as well
    mri_surf2surf --srcsubject fsaverage --trgsubject fsaverage6 --sval ${hemi}.somatomotor-mask.func.gii  --hemi lh --tval ${hemi}.somatomotor-mask.fsaverage6.func.gii 
    
done

freeview -f /usr/local/freesurfer/subjects/fsaverage/surf/rh.inflated:overlay=rh.somatomotor-mask.func.gii
freeview -f /usr/local/freesurfer/subjects/fsaverage/surf/lh.inflated:overlay=lh.somatomotor-mask.func.gii

LABELSDIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage6
mkdir $LABELSDIR
cd $LABELSDIR
mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi rh --annotation aparc.a2009s --subject fsaverage6
mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi lh --annotation aparc.a2009s --subject fsaverage6

# produce also individual labels in fsaverag6 space

for hemi in lh rh; do
 
  for label in $LABELS; do
    echo $label
    mris_convert --label ${hemi}.${label}.label ${label}\
        $SUBJECTS_DIR/fsaverage6/surf/${hemi}.white \
        ${hemi}.${label}.label.gii

    $HCPDIR/wb_command -gifti-label-to-roi ${hemi}.${label}.label.gii -name ${label} ${hemi}.${label}.func.gii
      

  done
done
#$HCPDIR/wb_command -metric-resample -to-nifti \
#    ${hemi}.mask.func.gii \
#    ${hemi}.mask.nii.gz


