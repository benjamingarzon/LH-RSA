#!/bin/sh

#labels of interest
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject

LABELSDIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage

VBM_TEMPLATE=/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/mask.nii.gz
fMRI_TEMPLATE=/home/benjamin.garzon/Data/LeftHand/Lund1/fmriprep/fmriprep/sub-lue1202/ses-1/func/sub-lue1202_ses-1_task-sequence_run-01_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz

mkdir $LABELSDIR
cd $LABELSDIR

# use Desikan
#mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi lh --annotation aparc.a2009s --subject fsaverage
#mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR  --hemi rh --annotation aparc.a2009s --subject fsaverage

#LABELS="S_front_inf G_postcentral S_central S_postcentral G_precentral S_front_sup Pole_occipital S_intrapariet_and_P_trans G_temporal_middle G_and_S_cingul-Mid-Ant S_intrapariet_and_P_trans G_front_sup S_precentral-sup-part G_parietal_sup G_and_S_paracentral S_front_sup S_cingul-Marginalis G_and_S_cingul-Mid-Ant G_and_S_cingul-Mid-Post G_precuneus G_front_middle S_precentral-inf-part"

# use Glasser
#mri_segstats --i fsaverage_aparc/rh.somatomotor-mask.func.gii --annot fsaverage rh ./GlasserParc/rh.HCP-MMP1.annot --sum parcels
#for i in `cat parcels | grep R_ | grep -v "0.0000     0.0000     0.0000     0.0000" | sed 's/\s\+/ /g' | cut -d 'R' -f2-5| cut -d' ' -f1`; do echo -n "R${i} "; done > rh.parcels; rm parcels
#mri_segstats --i fsaverage_aparc/rh.somatomotor-mask.func.gii --annot fsaverage lh ./GlasserParc/lh.HCP-MMP1.annot --sum parcels
#for i in `cat parcels | grep _ | grep -v "0.0000     0.0000     0.0000     0.0000" | sed 's/\s\+/ /g' | cut -d 'L' -f2-5| cut -d' ' -f1`; do echo -n "L${i} "; done > lh.parcels; rm parcels
#cat rh.parcels | sed 's/R_/L_/g' > lh.parcels 
cat ~/Software/LeftHand/masks/rh.parcels.txt | sed 's/R_/L_/g' > ~/Software/LeftHand/masks/lh.parcels.txt

mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi lh --annotation HCP-MMP1 --subject fsaverage
mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR  --hemi rh --annotation HCP-MMP1 --subject fsaverage

LABELS="`cat ~/Software/LeftHand/masks/rh.parcels.txt` `cat ~/Software/LeftHand/masks/lh.parcels.txt`"

for hemi in lh rh; do
 
  for label in $LABELS; do
    echo $label
    if [ ! -e ${hemi}.${label}.label ]; then
       continue
    fi
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
#LABELS="S_precentral-sup-part G_postcentral S_central G_precentral S_postcentral S_intrapariet_and_P_trans G_front_sup G_parietal_sup G_and_S_paracentral S_front_sup S_cingul-Marginalis G_and_S_cingul-Mid-Ant G_and_S_cingul-Mid-Post G_precuneus G_front_middle S_precentral-inf-part"

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

    # resample to MNI space as well     	
    mri_cor2label --i ./${hemi}.somatomotor-mask.func.gii \
    --id 1 --l ./${hemi}.somatomotor-mask.MNI.label \
    --surf fsaverage $hemi white 

   mri_label2label --srclabel ./${hemi}.somatomotor-mask.MNI.label \
   --srcsubject fsaverage \
   --trgsubject cvs_avg35_inMNI152 \
   --trglabel ./${hemi}.somatomotor-mask.MNI.label \
   --hemi ${hemi} --regmethod surface
    
   mri_label2vol --label ./${hemi}.somatomotor-mask.MNI.label \
    --regheader T1.mgz --subject cvs_avg35_inMNI152 --fill-ribbon --hemi $hemi \
    --o ./${hemi}.somatomotor-mask.MNI.nii.gz  

done

    fslmaths lh.somatomotor-mask.MNI.nii.gz -add rh.somatomotor-mask.MNI.nii.gz \
    -bin -kernel boxv 7 -dilF -bin somatomotor-mask.MNI.nii.gz

    # for VBM 
    mri_vol2vol --mov ./somatomotor-mask.MNI.nii.gz --targ $VBM_TEMPLATE --regheader --o ./somatomotor-mask.MNI.VBM.nii.gz --nearest
    fslmaths somatomotor-mask.MNI.VBM.nii.gz -mas $VBM_TEMPLATE somatomotor-mask.MNI.VBM.nii.gz

    # for MRI
    mri_vol2vol --mov ./somatomotor-mask.MNI.nii.gz --targ $VBM_TEMPLATE --regheader --o ./somatomotor-mask.MNI.fMRI.nii.gz --nearest
    fslmaths somatomotor-mask.MNI.fMRI.nii.gz -mas $fMRI_TEMPLATE somatomotor-mask.MNI.fMRI.nii.gz

freeview somatomotor-mask.MNI.fMRI.nii.gz -f /usr/local/freesurfer/subjects/fsaverage/surf/rh.inflated:overlay=rh.somatomotor-mask.func.gii:annot=../GlasserParc/rh.HCP-MMP1.annot
freeview somatomotor-mask.MNI.VBM.nii.gz -f /usr/local/freesurfer/subjects/fsaverage/surf/lh.inflated:overlay=lh.somatomotor-mask.func.gii:annot=../GlasserParc/lh.HCP-MMP1.annot

freeview -f /usr/local/freesurfer/subjects/fsaverage/surf/rh.inflated:overlay=rh.somatomotor-mask.func.gii -viewport 3D -cam azimuth 180 elevation 0 -zoom 1.3 -ss rh.mask_medial.png
freeview -f /usr/local/freesurfer/subjects/fsaverage/surf/lh.inflated:overlay=lh.somatomotor-mask.func.gii -viewport 3D -cam azimuth 180 elevation 0 -zoom 1.3 -ss lh.mask_lateral.png
freeview -f /usr/local/freesurfer/subjects/fsaverage/surf/rh.inflated:overlay=rh.somatomotor-mask.func.gii -viewport 3D -cam azimuth 0 elevation 0 -zoom 1.3 -ss rh.mask_lateral.png
freeview -f /usr/local/freesurfer/subjects/fsaverage/surf/lh.inflated:overlay=lh.somatomotor-mask.func.gii -viewport 3D -cam azimuth 0 elevation 0 -zoom 1.3 -ss lh.mask_medial.png

pngappend rh.mask_medial.png + lh.mask_medial.png medial.png
pngappend rh.mask_lateral.png + lh.mask_lateral.png lateral.png
pngappend medial.png - lateral.png mask.png

rm rh.mask_medial.png lh.mask_medial.png medial.png rh.mask_lateral.png lh.mask_lateral.png lateral.png
cp mask.png ~/Software/LeftHand/masks/

fslmaths /home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/rh.mask.nii.gz /home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/rh.mask.nii.gz -odt char
fslmaths /home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/lh.mask.nii.gz /home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/lh.mask.nii.gz -odt char
ln -s /home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/rh.mask.nii.gz /home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results/rh.mask.nii.gz
ln -s /home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/lh.mask.nii.gz /home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results/lh.mask.nii.gz


LABELSDIR=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage6
mkdir $LABELSDIR
cd $LABELSDIR
#mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi rh --annotation aparc.a2009s --subject fsaverage6
#mri_annotation2label --sd $SUBJECTS_DIR --outdir $LABELSDIR --hemi lh --annotation aparc.a2009s --subject fsaverage6

# produce also individual labels in fsaverage6 space

for hemi in lh rh; do
 
  for label in $LABELS; do
    echo $label
    if [ ! -e ../fsaverage/${hemi}.${label}.label ]; then
       continue
    fi

    mri_label2label --srclabel ../fsaverage/${hemi}.${label}.label \
    --srcsubject fsaverage \
    --trgsubject fsaverage6 \
    --trglabel ./${hemi}.${label}.label \
    --hemi ${hemi} --regmethod surface

    mris_convert --label ${hemi}.${label}.label ${label}\
        $SUBJECTS_DIR/fsaverage6/surf/${hemi}.white \
        ${hemi}.${label}.label.gii

    $HCPDIR/wb_command -gifti-label-to-roi ${hemi}.${label}.label.gii -name ${label} ${hemi}.${label}.func.gii
      

  done
done
#$HCPDIR/wb_command -metric-resample -to-nifti \
#    ${hemi}.mask.func.gii \
#    ${hemi}.mask.nii.gz


