#/usr/bin/sh
###############################################
# Convert to BIDS
###############################################

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lundpilot1
conda activate lhenv2

cd $HOMEDIR/

convert(){

SUBJECT=$1 
SESSION=$2
HEUDICONV_FILE=heudiconv_files/heuristic_conv-sub-${SUBJECT}-ses-${SESSION}.py
heudiconv -d '7T033/7T033{subject}_*/*/Dicom/*/*.dcm' -s ${SUBJECT} -f $HEUDICONV_FILE -o data_BIDS/sub-${SUBJECT}/ses-${SESSION} -b --overwrite
rm -r data_BIDS/sub-${SUBJECT}/ses-${SESSION}/sourcedata

}

convert 104 1 &
convert 105 1 &
convert 106 1 &
convert 105 2


