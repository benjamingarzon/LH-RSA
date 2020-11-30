#!/bin/sh

NPROCS=10
USENEWTEMPLATE=0

VBM_DIR="/home/benjamin.garzon/Data/LeftHand/Lund1/vbm"
BIDS_DIR="/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS"
CAT_DIR="/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/data"
CAT_DIR="/home/benjamin.garzon/Data/LeftHand/Lund1/spm12/data"
STRUC_DIR="$VBM_DIR/struc"
BOUNDINGBOX="0 257 0 320 40 280"
mkdir $STRUC_DIR
scan_list_file="$VBM_DIR/scans.txt"
subject_list_file="$VBM_DIR/subjects.txt"

cd $VBM_DIR
printf '%s\n' sub-*/ses* > $scan_list_file
printf '%s\n' sub-* > $subject_list_file

myimage=T1wtotemplate

scanlist=""
scans=`cat $scan_list_file`
#create links

# for CAT12
rm boundingboxes.txt
for scan in $scans; do
    NAME=`echo $scan | sed 's@/ses-@.@g'`
    echo $NAME
    rm $CAT_DIR/${NAME}_MP2RAGE.nii
    fslmaths $VBM_DIR/${scan}/brainmasknative.nii.gz -dilF -kernel boxv 5 $CAT_DIR/mask.nii.gz 
    BOUNDINGBOX=`fslstats  $CAT_DIR/mask.nii.gz -w`
    echo $BOUNDINGBOX >> boundingboxes.txt
done

for scan in $scans; do
    NAME=`echo $scan | sed 's@/ses-@.@g'`
    echo $NAME
    rm $CAT_DIR/${NAME}_MP2RAGE.nii
    fslmaths $VBM_DIR/${scan}/brainmasknative.nii.gz -dilF -kernel boxv 5 $CAT_DIR/mask.nii.gz 
    BOUNDINGBOX=`fslstats  $CAT_DIR/mask.nii.gz -w`
    rm $CAT_DIR/mask.nii.gz
    fslroi $BIDS_DIR/${scan}/anat/MP2RAGEpos.nii.gz $CAT_DIR/${NAME}_MP2RAGE.nii.gz $BOUNDINGBOX
    fslreorient2std $CAT_DIR/${NAME}_MP2RAGE.nii.gz $CAT_DIR/${NAME}_MP2RAGE.nii.gz
    gunzip $CAT_DIR/${NAME}_MP2RAGE.nii.gz
done

