#/usr/bin/bash
###############################################
# Convert to BIDS
###############################################

HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
PROGDIR=/home/benjamin.garzon/Software/LeftHand/process/
SCANLIST=$HOMEDIR/dicoms/Scan_list_wave1_complete.csv
SCANLIST=$HOMEDIR/dicoms/Scan_list_wave1_missing.csv
export HEUDICONV_FILE=/home/benjamin.garzon/Data/LeftHand/Lund1/heudiconv_files/heuristic_conv.py

#clean up .heudiconv 
rm -r $HOMEDIR/data_BIDS/.heudiconv

conda init /usr/bin/bash 
conda activate lhenv2

cd $HOMEDIR/

convert(){

SUBJECT=$1 
SESSION=$2
SUFFIX=$3

if  [ ! -e "dicoms/$SUBJECT/${SUFFIX}" ] && [ -e "dicoms/$SUBJECT/${SUFFIX}.tar.gz" ]; then
    echo "Uncompressing  $SUBJECT $SESSION"
    tar -xzf dicoms/$SUBJECT/${SUFFIX}.tar.gz
fi
# take the chance to clean up
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*dummy*
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*SmartBrain*
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*Aligned*
rm -r 7T033{subject}_${SUFFIX}/*/Dicom/*Survey*

echo heudiconv -d "dicoms/{subject}/${SUFFIX}/Dicom/*/*.dcm" -s ${SUBJECT} -ss ${SESSION} -f $HEUDICONV_FILE -o data_BIDS/ -b --overwrite
heudiconv -d "dicoms/{subject}/${SUFFIX}/Dicom/*/*.dcm" -s ${SUBJECT} -ss ${SESSION} -f $HEUDICONV_FILE -o data_BIDS/ -b --overwrite

# compress and clean when finished
echo "Compressing $SUBJECT $SESSION"
tar -czf dicoms/$SUBJECT/${SUFFIX}.tar.gz dicoms/$SUBJECT/${SUFFIX} 
rm -r dicoms/$SUBJECT/${SUFFIX}
rm -r data_BIDS/sourcedata/sub-${SUBJECT}/ses-${SESSION}
}

rm $HOMEDIR/dicoms/missing.csv

tail -n +2 $SCANLIST | while read line
do
SUBJECT=`echo $line | cut -f1 -d','`
SESSION=`echo $line | cut -f2 -d','`
DIR=`echo $line | cut -f3 -d','`
TAB='\t' 
head -n 1 $SCANLIST | sed "s/,/${TAB}/g" > $HOMEDIR/dicoms/curr_subject.csv
echo $line | sed "s/,/${TAB}/g"  >> $HOMEDIR/dicoms/curr_subject.csv

    if [ ! -e "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/anat" ] || [ ! -e "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/fmap" ] || [ ! -e "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/func" ]; then
        echo Converting $SUBJECT $SESSION $DIR
        convert $SUBJECT $SESSION $DIR
    fi

    # keep track of missing data
   if [ ! -e "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/anat" ]; then
        echo Anatomical missing: "sub-$SUBJECT ses-$SESSION" >> $HOMEDIR/dicoms/missing.csv
     else 
     if [ `ls "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/anat" | wc -l` -lt 8 ]; then 
        echo Anatomical missing: "sub-$SUBJECT ses-$SESSION" >> $HOMEDIR/dicoms/missing.csv
     fi
   fi

   if [ ! -e "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/fmap" ]; then
        echo Fieldmap missing: "sub-$SUBJECT ses-$SESSION" >> $HOMEDIR/dicoms/missing.csv
     else 
     if [ `ls "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/fmap" | wc -l` -lt 4 ]; then 
        echo Fieldmap missing: "sub-$SUBJECT ses-$SESSION" >> $HOMEDIR/dicoms/missing.csv
     fi
   fi

   if [ ! -e "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/func" ]; then
        echo Func runs missing: "sub-$SUBJECT ses-$SESSION" >> $HOMEDIR/dicoms/missing.csv
     else 
     if [ `ls "$HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/func" | wc -l` -lt 15 ]; then 
        echo Func runs missing: "sub-$SUBJECT ses-$SESSION. Only found `ls $HOMEDIR/data_BIDS/sub-$SUBJECT/ses-$SESSION/func/*.nii.gz | wc -l`" >> $HOMEDIR/dicoms/missing.csv
     fi
   fi


rm $HOMEDIR/dicoms/curr_subject.csv
done


#convert lue001 1 20190426 &
#convert lue001 2 20190506 &
#convert lue001 3 20190513 &