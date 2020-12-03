#!/bin/sh
HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

WD=$1
METRIC=$2
DISTANCE=$3
RADIUS=$4
OUTPUT=$5
THR=$6

EXTREMA=$WD/extrema.nii.gz

#if [ ! -e $EXTREMA ]; then
$HCPDIR/wb_command -volume-extrema \
      $METRIC \
      $DISTANCE \
      $EXTREMA -only-maxima -threshold 0 $THR
#fi

if [ `$HCPDIR/wb_command -volume-stats $EXTREMA -reduce MAX` -eq '0' ]; then 
echo "No clusters!"
cp $EXTREMA ${OUTPUT}.nii.gz
cp $EXTREMA ${OUTPUT}.all.nii.gz
exit 1
fi

#if [ ! -e ${OUTPUT}.nii.gz ]; then
$HCPDIR/wb_command -volume-rois-from-extrema \
      $EXTREMA \
      $RADIUS \
      ${OUTPUT}.nii.gz
#fi

#cp $EXTREMA ${OUTPUT}.nii.gz

$HCPDIR/wb_command -volume-reduce \
      ${OUTPUT}.nii.gz \
      INDEXMAX \
      ${OUTPUT}.aux.nii.gz

$HCPDIR/wb_command -volume-math \
      '(x>1)*x + y' \
      ${OUTPUT}.all.nii.gz \
      -var x ${OUTPUT}.aux.nii.gz \
      -var y ${OUTPUT}.nii.gz -subvolume 1

rm ${OUTPUT}.aux.nii.gz

$HCPDIR/wb_command -volume-stats \
      ${OUTPUT}.all.nii.gz\
      -reduce MAX
      

      