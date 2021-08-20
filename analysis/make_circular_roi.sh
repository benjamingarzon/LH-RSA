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
MINCLUSSIZE=100.0

EXTREMA=$WD/extrema.func.gii
#if [ ! -e $EXTREMA ]; then

      
$HCPDIR/wb_command -metric-extrema \
      $SURFACE \
      $METRIC \
      $DISTANCE \
      $EXTREMA -only-maxima -threshold 0 $THR # -consolidate-mode
#fi

$HCPDIR/wb_command -metric-extrema \
      $SURFACE \
      $METRIC \
      $DISTANCE \
      $EXTREMA -only-maxima -threshold 0 $THR # -consolidate-mode

if [ `$HCPDIR/wb_command -metric-stats $EXTREMA -reduce MAX` -eq '0' ]; then 
echo "No clusters!"
cp $EXTREMA ${OUTPUT}.func.gii
cp $EXTREMA ${OUTPUT}.all.func.gii
$HCPDIR/wb_command -metric-convert -to-nifti $EXTREMA ${OUTPUT}.nii.gz
exit 1
fi

#if [ ! -e ${OUTPUT}.func.gii ]; then
$HCPDIR/wb_command -metric-rois-from-extrema \
      $SURFACE \
      $EXTREMA \
      $RADIUS \
      ${OUTPUT}.func.gii 
#fi


$HCPDIR/wb_command -metric-reduce \
      ${OUTPUT}.func.gii \
      INDEXMAX \
      ${OUTPUT}.aux.func.gii 

$HCPDIR/wb_command -metric-math \
      '(x>1)*x + y' \
      ${OUTPUT}.all.func.gii \
      -var x ${OUTPUT}.aux.func.gii \
      -var y ${OUTPUT}.func.gii -column 1

rm ${OUTPUT}.aux.func.gii
 
NCLUS=`$HCPDIR/wb_command -metric-stats \
      ${OUTPUT}.all.func.gii \
      -reduce MAX`
echo $NCLUS
#for i in `seq $NCLUS`; do      
#   CLUSINDEX=`$HCPDIR/wb_command -metric-stats \
#         $EXTREMA \
#         -reduce INDEXMAX -column $i`
#         echo $i $CLUSINDEX
#done

$HCPDIR/wb_command -metric-convert -to-nifti ${OUTPUT}.func.gii ${OUTPUT}.nii.gz

mri_convert ${OUTPUT}.func.gii ${OUTPUT}.func.gii
mri_convert ${OUTPUT}.all.func.gii ${OUTPUT}.all.func.gii
rm $EXTREMA

mri_surfcluster --in $METRIC --thmin $THR --minarea $MINCLUSSIZE --hemi $HEMI \
--sum ${OUTPUT}.clusters.sum --subject fsaverage6 --ocn ${OUTPUT}.clusters.func.gii
cat ${OUTPUT}.clusters.sum
#system('if [ `cat clusters.sum | wc -l` -eq 34 ]; then rm clusters.sum; fi') # no clusters
  