#!/usr/bin/sh

DEPTHS="0.20 0.40 0.60 0.80" 
for i in $DEPTHS; do 
  for j in 1 2 3 4 5; do 
     echo --------------------------------------------------------
     echo Thickness wave $j depth $i
     echo --------------------------------------------------------
    ./prepare_fs_analysis.sh T1 $i $j &
  done 
done 

exit 1 
# run:
for j in 1 2 3 4 5; do 
echo --------------------------------------------------------
echo Thickness wave $j
echo --------------------------------------------------------
./prepare_fs_analysis.sh thickness 0 $j &
done
exit 1

