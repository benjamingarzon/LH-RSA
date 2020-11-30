#!/usr/bin/sh


DEPTHS="-0.20 0.00 0.20 0.40 0.60 0.80" 
for i in $DEPTHS; do 
  echo --------------------------------------------------------
  echo T1 depth $i
  echo --------------------------------------------------------

  seq 5 | xargs -P 5 -i ./prepare_fs_analysis.sh T1 $i {}
#  seq 5 | xargs -i ./prepare_fs_analysis.sh T1 $i {}


done 

exit 1 

# run:
# seq 5 | xargs -P 5 -i ./prepare_fs_analysis.sh thickness 0 {}

exit 1

