#/usr/bin/bash
###############################################
# Convert to BIDS
###############################################

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lundpilot1
PROGDIR=/home/benjamin.garzon/Software/LeftHand/process/

conda init /usr/bin/bash 
conda activate lhenv2

cd $HOMEDIR/

convert(){

SUBJECT=$1 
SESSION=$2
SUFFIX=$3
#python $PROGDIR/fix_dicoms.py 7T033/7T033${SUBJECT}
HEUDICONV_FILE=heudiconv_files/heuristic_conv-sub-${SUBJECT}-ses-${SESSION}.py

heudiconv -d "7T033/7T033{subject}_${SUFFIX}/*/Dicom/*/*.dcm" -s ${SUBJECT} -ss ${SESSION} -f $HEUDICONV_FILE -o data_BIDS/ -b --overwrite
rm -r data_BIDS/sourcedata/sub-${SUBJECT}/ses-${SESSION}

}

convert 105 3 20190404 &
#convert 107 1 20190326 &
#sleep 10
#convert 107 2 20190329 &
#sleep 10
#convert 103 3 20190401

#convert 106 1 20190322 

#convert 103 1 20190225 &
#sleep 10
#convert 103 2 20190312 &
#sleep 10
#convert 104 1 20190228 &
#sleep 10
#convert 105 1 20190228 &
#sleep 10
#convert 105 2 20190311 
