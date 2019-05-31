#/usr/bin/bash
###############################################
# Convert to BIDS
###############################################

#HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lundpilot1
HOMEDIR=/home/share/MotorSkill
PROGDIR=/home/benjamin.garzon/Software/LeftHand/process/
SCANLIST=$HOMEDIR/dicoms/Scan_list_wave1x.csv

cp $HOMEDIR/dicoms/Scan_list_wave1.csv $HOMEDIR/dicoms/Scan_list.csv

conda init /usr/bin/bash 
conda activate lhenv2

cd $HOMEDIR/

convert(){

SUBJECT=$1 
SESSION=$2
SUFFIX=$3
#python $PROGDIR/fix_dicoms.py 7T033/7T033${SUBJECT}
#HEUDICONV_FILE=heudiconv_files/heuristic_conv-sub-${SUBJECT}-ses-${SESSION}.py
HEUDICONV_FILE=heudiconv_files/heuristic_conv.py
if  [ ! -e "dicoms/$SUBJECT/${SUFFIX}" ] && [ -e "dicoms/$SUBJECT/${SUFFIX}.tar.gz" ]; then
    tar -xvzf dicoms/$SUBJECT/${SUFFIX}.tar.gz
fi
# take the chance to clean up
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*dummy*
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*SmartBrain*
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*Aligned*
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*Survey*

heudiconv -d "dicoms/{subject}/${SUFFIX}/*/Dicom/*/*.dcm" -s ${SUBJECT} -ss ${SESSION} -f $HEUDICONV_FILE -o data_BIDS/ -b --overwrite
#heudiconv -d "dicoms/{subject}/7T033???_${SUFFIX}/*/Dicom/*/*.dcm" -s ${SUBJECT} -ss ${SESSION} -f $HEUDICONV_FILE -o data_BIDS/ -b --overwrite
# compress and clean when finished
tar -cvzf dicoms/$SUBJECT/${SUFFIX}.tar.gz dicoms/$SUBJECT/${SUFFIX} 
rm -r dicoms/$SUBJECT/${SUFFIX}
rm -r data_BIDS/sourcedata/sub-${SUBJECT}/ses-${SESSION}
}

tail -n +2 $SCANLIST | while read line
do
SUBJECT=`echo $line | cut -f1 -d' '`
SESSION=`echo $line | cut -f2 -d' '`
DIR=`echo $line | cut -f3 -d' '`

    if [ ! -e "data_BIDS/sub-$SUBJECT/ses-$SESSION" ]; then
        echo Converting $SUBJECT $SESSION $DIR 
        convert $SUBJECT $SESSION $DIR
    fi
done


#convert lue001 1 20190426 &
#convert lue001 2 20190506 &
#convert lue001 3 20190513 &

#convert 109 1 20190426 &
#convert 109 2 20190506 &
#convert 109 3 20190513 &

#convert 105 3 20190404 &
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
