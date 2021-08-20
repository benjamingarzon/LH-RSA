#!/bin/sh
# show surface map

#SURF_L=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/lh.inflated
#SURF_R=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/rh.inflated
#CURV_L=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/lh.curv
#CURV_R=/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/rh.curv

WD=$1
OVERLAY_L=$2
OVERLAY_R=$3
THR_L=$4
THR_H=$5
FIG=$6
SURF_L=$7
SURF_R=$8
ANNOT_L=$9
ANNOT_R=${10}

cd $WD

freeview -f $SURF_L:overlay=$OVERLAY_L:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse':overlay_opacity=0.6:$ANNOT_L \
-viewport 3D -cam azimuth 180 elevation 0 -zoom 1.7 -ss ${FIG}_L_medial 

freeview -f $SURF_L:overlay=$OVERLAY_L:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse':overlay_opacity=0.6:$ANNOT_L \
-viewport 3D -cam azimuth 0 elevation 0 -zoom 1.5 -ss ${FIG}_L_lateral

freeview -f $SURF_R:overlay=$OVERLAY_R:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse':overlay_opacity=0.6:$ANNOT_R \
-viewport 3D -colorscale -cam azimuth 180 elevation 0 -zoom 1.5 -ss ${FIG}_R_lateral

freeview -f $SURF_R:overlay=$OVERLAY_R:overlay_threshold=$THR_L,$THR_H:overlay_color='colorwheel','inverse':overlay_opacity=0.6:$ANNOT_R \
-viewport 3D -cam azimuth 0 elevation 0 -zoom 1.7 -ss ${FIG}_R_medial

pngappend ${FIG}_L_medial.png + ${FIG}_R_medial.png ${FIG}_medial.png
pngappend ${FIG}_L_lateral.png + ${FIG}_R_lateral.png ${FIG}_lateral.png
pngappend ${FIG}_medial.png - ${FIG}_lateral.png ${FIG}.png

rm ${FIG}_medial.png ${FIG}_lateral.png

if [ -e  ${FIG}.png ]; then 
rm ${FIG}_L_medial.png ${FIG}_R_medial.png ${FIG}_L_lateral.png ${FIG}_R_lateral.png
fi
