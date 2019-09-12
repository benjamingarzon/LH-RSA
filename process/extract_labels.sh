#!/bin/sh

mri_annotation2label --sd /usr/local/freesurfer/subjects/ --outdir /home/benjamin.garzon/Data/LeftHand/Lund1/labels --hemi lh --annotation aparc.a2009s --subject fsaverage
mri_annotation2label --sd /usr/local/freesurfer/subjects/ --outdir /home/benjamin.garzon/Data/LeftHand/Lund1/labels --hemi rh --annotation aparc.a2009s --subject fsaverage

cd /home/benjamin.garzon/Data/LeftHand/Lund1/labels
#labels of interest
LABELS="S_front_inf S_postcentral G_precentral S_front_sup Pole_occipital S_intrapariet_and_P_trans G_temporal_middle G_and_S_cingul-Mid-Ant"

HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/
for label in $LABELS; do
    echo $label
    mris_convert --label lh.${label}.label ${label}\
        /usr/local/freesurfer/subjects/fsaverage/surf/lh.white \
        lh.${label}.label.gii
    $HCPDIR/wb_command -gifti-label-to-roi lh.${label}.label.gii -name ${label} lh.${label}.func.gii
        
    mris_convert --label rh.${label}.label ${label}\
        /usr/local/freesurfer/subjects/fsaverage/surf/rh.white \
        rh.${label}.label.gii
    $HCPDIR/wb_command -gifti-label-to-roi rh.${label}.label.gii -name ${label} rh.${label}.func.gii

done


