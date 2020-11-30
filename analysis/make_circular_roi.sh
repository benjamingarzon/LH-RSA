#!/bin/sh
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

WD=$1
SURFACE=$2
METRIC=$3
DISTANCE=$4
RADIUS=$5
OUTPUT=$6
THR=$7

EXTREMA=$WD/extrema.func.gii
#if [ ! -e $EXTREMA ]; then
$HCPDIR/wb_command -metric-extrema \
      $SURFACE \
      $METRIC \
      $DISTANCE \
      $EXTREMA -only-maxima -threshold 0 $THR
#fi

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

$HCPDIR/wb_command -metric-stats \
      ${OUTPUT}.all.func.gii \
      -reduce MAX
      
$HCPDIR/wb_command -metric-convert -to-nifti ${OUTPUT}.func.gii ${OUTPUT}.nii.gz

rm $EXTREMA