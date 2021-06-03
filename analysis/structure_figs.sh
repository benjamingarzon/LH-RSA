#!/bin/sh
# create figs for structural analyses

HOMEDIR=/home/benjamin.garzon
FIGS_DIR=$HOMEDIR/Data/LeftHand/Lund1/figs/structure
RH_INFLATED=$HOMEDIR/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/rh.inflated
LH_INFLATED=$HOMEDIR/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/lh.inflated
RH_LABEL=$HOMEDIR/Data/LeftHand/Lund1/labels/fsaverage/rh.somatomotor-mask.label
LH_LABEL=$HOMEDIR/Data/LeftHand/Lund1/labels/fsaverage/lh.somatomotor-mask.label

#mri_cor2label --i $HOMEDIR/Data/LeftHand/Lund1/labels/fsaverage/rh.somatomotor-mask.func.gii --id 1 --l $HOMEDIR/Data/LeftHand/Lund1/labels/fsaverage/rh.somatomotor-mask.label --sd $HOMEDIR/Data/LeftHand/Lund1/labels/subject/ --surf fsaverage rh white
#mri_cor2label --i $HOMEDIR/Data/LeftHand/Lund1/labels/fsaverage/lh.somatomotor-mask.func.gii --id 1 --l $HOMEDIR/Data/LeftHand/Lund1/labels/fsaverage/lh.somatomotor-mask.label --sd $HOMEDIR/Data/LeftHand/Lund1/labels/subject/ --surf fsaverage lh white

# 1. Reliability cortical thickness
WD=$HOMEDIR/Data/LeftHand/Lund1/freesurfer/results/tests/thickness
$HOMEDIR/Software/LeftHand/analysis/show_surface.sh $FIGS_DIR $WD/reliability.lh/ICC.func.gii $WD/reliability.rh/ICC.func.gii 0.5 1 reliability-thickness $LH_INFLATED $RH_INFLATED $LH_LABEL $RH_LABEL

WD=$HOMEDIR/Data/LeftHand/Lund1/freesurfer/results/tests/thickness
$HOMEDIR/Software/LeftHand/analysis/show_surface.sh $FIGS_DIR $WD/reliability-intervention.lh/ICC.func.gii $WD/reliability-intervention.rh/ICC.func.gii 0.5 1 reliability-thickness-intervention $LH_INFLATED $RH_INFLATED $LH_LABEL $RH_LABEL

# 2. Reliability VBM
cd /home/benjamin.garzon/Data/LeftHand/Lund1/cat12/tests/reliability
fsleyes --scene ortho --worldLoc -0.37494659423828125 -17.62505340576172 18.37494659423828 --displaySpace /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz --xcentre  0.00000  0.00000 --ycentre  0.00000  0.00000 --zcentre  0.00000  0.00000 --xzoom 100.0 --yzoom 100.0 --zzoom 100.0 --layout horizontal --hideCursor --bgColour 0.0 0.0 0.0 --fgColour 1.0 1.0 1.0 --cursorColour 0.0 1.0 0.0 --showColourBar --colourBarLocation right --colourBarLabelSide top-left --colourBarSize 100.0 --labelSize 12 --performance 3 --movieSync /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz --name "MNI152_T1_1mm_brain" --overlayType volume --alpha 100.0 --brightness 50.0 --contrast 50.0 --cmap greyscale --negativeCmap greyscale --displayRange 0.0 8364.0 --clippingRange 0.0 8447.64 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 100 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 /home/benjamin.garzon/Data/LeftHand/Lund1/cat12/tests/reliability/ICC.nii.gz --name "ICC" --overlayType volume --alpha 55.33333333566164 --brightness 37.4427772923199 --contrast 74.96185152821327 --cmap brain_colours_x_rain --negativeCmap greyscale --displayRange 0.5 1.0 --clippingRange 0.5 1.008461149930954 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 100 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 /home/benjamin.garzon/Data/LeftHand/Lund1/cat12/mask.nii.gz --name "somatomotor-mask.MNI" --overlayType label --alpha 100.0 --brightness 50.0 --contrast 50.0 --lut mgh-cma-freesurfer --outline --outlineWidth 3 --volume 0# add histograms of reliability



#3. Univariate thickness
# run roi_analyses.R

# 4. Univariate VBM
# roi_analyses.R

# 5. Multivariate thickness
WD=$HOMEDIR/Data/LeftHand/Lund1/freesurfer/results_unsmoothed/tests/multivariate-tess/
$HOMEDIR/Software/LeftHand/analysis/show_surface.sh $FIGS_DIR/multivar-thickness $WD/lh.accuracy.map $WD/rh.accuracy.map 0.55 0.65 multivar-tess-thickness $LH_INFLATED $RH_INFLATED $LH_LABEL $RH_LABEL

# 6. Coefficients
# 7. Violin plots for accuracy

#v8. Multivariate VBM
WD=$HOMEDIR/Data/LeftHand/Lund1/cat12//tests/multivariate-tess/
$HOMEDIR/Software/LeftHand/analysis/show_surface.sh $FIGS_DIR/multivar-VBM $WD/lh.accuracy.map $WD/rh.accuracy.map 0.55 0.65 multivar-tess-vbm $LH_INFLATED $RH_INFLATED $LH_LABEL $RH_LABEL
# 9. Coefficients
# 10. Violin plots for accuracy
