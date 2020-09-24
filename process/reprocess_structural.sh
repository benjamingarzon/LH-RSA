#!/bin/sh
# reprocessing structural data
# redo those subjects that had 

#
./run_structural_pipeline.sh lue3205 "2 3 4 5 6"

exit 1
./run_structural_pipeline.sh lue3105 "2 3 4 5 6 7"
./run_structural_pipeline.sh lue2207 "1"
./run_structural_pipeline.sh lue2102 "1 2 3 4 5 7" 
./run_structural_pipeline.sh lue2103 "1 2 3 4 6 7" 
./run_structural_pipeline.sh lue2104 "1 2 4 5 7"
./run_structural_pipeline.sh lue3102 "1 2 3 4 6 7"
./run_structural_pipeline.sh lue3104 "2 3 4 5 6 7"
./run_structural_pipeline.sh lue3106 "1"
./run_structural_pipeline.sh lue3203 "2 3"
./run_structural_pipeline.sh lue3205 "1 2 3 4 5 6"
./run_structural_pipeline.sh lue4201 "1 2 3 4 5 6"

exit 1

#./run_structural_pipeline.sh lue3102 "1 2 3 4 5 6 7"

#2102 6
#2103 5
#2104 3 6
#3102 5
#3104 1
#3105 1
#3202 ???

#3203 1
#3205 1 7 
#3207 --
#4102 --
#4103 --
#4207 --
#5106
#4201


# 1201.3, 1201.4, 1201.5 bad quality scan, not used for constructing base

#Visualize
SUBJECTS="sub-lue1101 sub-lue1103 sub-lue1104 sub-lue1105 sub-lue1106 sub-lue1107 sub-lue1201 sub-lue1202 sub-lue1203 sub-lue1204 sub-lue1205 sub-lue1206 sub-lue1207 sub-lue2101 sub-lue2102 sub-lue2103 sub-lue2104 sub-lue2105 sub-lue2106 sub-lue2107 sub-lue2201 sub-lue2202 sub-lue2203 sub-lue2205 sub-lue2206 sub-lue3101 sub-lue3102 sub-lue3103 sub-lue3104 sub-lue3107 sub-lue3201 sub-lue3202 sub-lue3203 sub-lue3204 sub-lue3205 sub-lue3206 sub-lue3207 sub-lue4101 sub-lue4102 sub-lue4103 sub-lue4104 sub-lue4105 sub-lue4106 sub-lue4107 sub-lue4201 sub-lue4202 sub-lue4203 sub-lue4204 sub-lue4205 sub-lue4206 sub-lue4207 sub-lue5101 sub-lue5102 sub-lue5103 sub-lue5104 sub-lue5105 sub-lue5106 sub-lue5107 sub-lue5201 sub-lue5202 sub-lue5203 sub-lue5204 sub-lue5205 sub-lue5206 sub-lue5207"

SUBJECTS="sub-lue2102 sub-lue2103 sub-lue2104 sub-lue3102 sub-lue3104 sub-lue3203 sub-lue3105 sub-lue3106 sub-lue2207" #sub-lue3205 sub-lue3206 sub-lue3207 sub-lue4102 sub-lue4103 sub-lue4207 sub-lue5106"
SUBJECTS="sub-lue1204 sub-lue3205"

for subject in $SUBJECTS; do 
  echo $subject; 
  #check surfaces
  freeview ${subject}.base/mri/T1.mgz -f ${subject}*base/surf/rh.pial $subject*base/surf/lh.pial $subject*base/surf/rh.white $subject*base/surf/lh.white 
#-viewport 'coronal' -ss pix/${subject}
done


