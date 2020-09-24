#!/bin/sh

#use FREESURFER 7 
#SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer
#export RECONALL=~/Software/LeftHand/process/recon-allT2
#VBM_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/vbm

#use FREESURFER 6 
if [ ]; then 
export FREESURFER_HOME=/usr/local/freesurfer6
export FSFAST_HOME=/usr/local/freesurfer6/fsfast
export MNI_DIR=/usr/local/freesurfer6/mni
#export RECONALL=/home/benjamin.garzon/Software/LeftHand/process/recon-allT2-6.0
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer6
VBM_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/vbm6
bash /usr/local/freesurfer6/SetUpFreeSurfer.sh
export PATH=`echo $PATH | sed 's%freesurfer%freesurfer6%g'` 
export FS_LICENSE=/usr/local/freesurfer/license.txt
else

export FREESURFER_HOME=/home/benjamin.garzon/Software/freesurfer
export FSFAST_HOME=/home/benjamin.garzon/Software/freesurfer/fsfast
export MNI_DIR=/home/benjamin.garzon/Software/freesurfer/mni
#export RECONALL=/home/benjamin.garzon/Software/LeftHand/process/recon-allT2-6.0
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer6
VBM_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/vbm6
bash $FREESURFER_HOME/SetUpFreeSurfer.sh
export PATH=`echo $PATH | sed 's%/usr/local/freesurfer%/home/benjamin.garzon/Software/freesurfer%g'` 
export FS_LICENSE=/home/benjamin.garzon/Software/freesurfer/license.txt
fi

which recon-all
# define vars
HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
EXPERT_FILE=/home/benjamin.garzon/Software/LeftHand/process/expert.opts

# remove myelin inhomogeneity

MAXPROCS=75
 
#SUBJECTS="lue1101 lue1102 lue1103 lue1104 lue1105 lue1106 lue1107 lue1201 lue1202 lue1203 lue1204 lue1205 lue1206 lue1207"
#SUBJECTS="lue3101 lue3102 lue3103 lue3104 lue3105 lue3106 lue3107 lue3201 lue3202 lue3203 lue3204 lue3205 lue3206 lue3207"
#SUBJECTS="lue2101 lue2102 lue2103 lue2104 lue2105 lue2106 lue2107 lue2201 lue2202 lue2203 lue2204 lue2205 lue2206 lue2207"
#SUBJECTS="lue4101 lue4102 lue4103 lue4104 lue4105 lue4106 lue4107 lue4201 lue4202 lue4203 lue4204 lue4205 lue4206 lue4207"
#SUBJECTS="lue5101 lue5102 lue5103 lue5104 lue5105 lue5106 lue5107 lue5201 lue5202 lue5203 lue5204 lue5205 lue5206 lue5207"
echo $1 $2
#SUBJECTS="lue4201"
if [ "$1" ]; then
SUBJECTS="$1"
fi

if [ "$2" ]; then
SESSIONS="$2"
else 
SESSIONS=""
fi

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS

OVERWRITE=1
DOFS=1
for SUBJECT in $SUBJECTS; do
#    if [ `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/*/anat/sub-${SUBJECT}_ses-*_MP2RAGE.nii.gz | wc -l` -eq $NSESSIONS ]; then    
       nice ./run_fs.sh $SUBJECTS_DIR $WD $VBM_DIR $SUBJECT $EXPERT_FILE $OVERWRITE $DOFS "$SESSIONS" > $HOMEDIR/logs/${SUBJECT}.fs.log &
#    fi
    
    sleep 10
    NPROCS=`pgrep run_fs |wc -w`
    while [ $NPROCS -gt $MAXPROCS ]; do
      echo "$NPROCS processes running"
      sleep 1000
      NPROCS=`pgrep run_fs |wc -w`
    done
done

#WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
#rm $WD/*/ses*/anat/*template*  $WD/*/ses*/anat/*weights* $WD/*/ses*/anat/*.lta $WD/*/ses*/anat/*native* $WD/*/ses*/anat/*restore* $WD/*/ses*/anat/*masked*

# rm -r $SUBJECTS_DIR/*.template


