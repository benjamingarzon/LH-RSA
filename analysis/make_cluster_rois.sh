#!/bin/sh
#HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/
HCPDIR=/home/xgarzb@GU.GU.SE/Software/workbench/bin_rh_linux64/

WD=$1
SURFACE=$2
METRIC=$3
DISTANCE=$4
RADIUS=$5
OUTPUT=$6
THR=$7
HEMI=$8
MINCLUSSIZE=50.0

EXTREMA=$WD/extrema.func.gii

mri_surfcluster --in $METRIC --thmin $THR --minarea $MINCLUSSIZE --hemi $HEMI \
--sum ${OUTPUT}.sum --subject fsaverage6 --ocn ${OUTPUT}.all.func.gii

mri_convert ${OUTPUT}.all.func.gii ${OUTPUT}.all.func.gii

$HCPDIR/wb_command -metric-convert -to-nifti ${OUTPUT}.all.func.gii ${OUTPUT}.nii.gz

cat ${OUTPUT}.sum
#if [ `cat clusters.sum | wc -l` -eq 34 ]; then rm clusters.sum; fi # no clusters
if [ `$HCPDIR/wb_command -metric-stats ${OUTPUT}.all.func.gii -reduce MAX` -eq '0' ]; then 
echo "No clusters!"
#exit 1
fi


