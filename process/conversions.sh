#/usr/bin/bash
###############################################
# Convert to BIDS
###############################################

export HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
PROGDIR=/home/benjamin.garzon/Software/LeftHand/process/

SCANLIST=$HOMEDIR/dicoms/ScanList_wave5.csv
export HEUDICONV_FILE=/home/benjamin.garzon/Data/LeftHand/Lund1/heudiconv_files/heuristic_conv_reim.py

#clean up .heudiconv 
rm -r $HOMEDIR/data_BIDS/.heudiconv

conda init /usr/bin/bash 
conda activate lhenv2

cd $HOMEDIR/

convert(){

cd $HOMEDIR/

SUBJECT=$1 
SESSION=$2
SUFFIX=$3

if  [ ! -e "$HOMEDIR/dicoms/$SUBJECT/${SUFFIX}" ] && [ -e "$HOMEDIR/dicoms/$SUBJECT/${SUFFIX}.tar.gz" ]; then
    echo "Uncompressing  $SUBJECT $SESSION $SUFFIX"
    #cd $HOMEDIR/dicoms/$SUBJECT/
    tar -xzf $HOMEDIR/dicoms/$SUBJECT/${SUFFIX}.tar.gz --directory $HOMEDIR/dicoms/$SUBJECT/
fi
# take the chance to clean up
rm -r $HOMEDIR/dicoms/$SUBJECT/${SUFFIX}/Dicom/*dummy*
rm -r $HOMEDIR/dicoms/$SUBJECT/${SUFFIX}/Dicom/*SmartBrain*
rm -r $HOMEDIR/dicoms/$SUBJECT/${SUFFIX}/Dicom/*Aligned*
rm -r $HOMEDIR/dicoms/$SUBJECT/${SUFFIX}/Dicom/*Survey*

echo heudiconv -d "$HOMEDIR/dicoms/{subject}/${SUFFIX}/Dicom/*/*.dcm" -s ${SUBJECT} -ss ${SESSION} -f $HEUDICONV_FILE -o $HOMEDIR/data_BIDS/ -b --overwrite
heudiconv -d "$HOMEDIR/dicoms/{subject}/${SUFFIX}/Dicom/*/*.dcm" -s ${SUBJECT} -ss ${SESSION} -f $HEUDICONV_FILE -o $HOMEDIR/data_BIDS/ -b --overwrite

# compress and clean when finished
echo "Compressing $SUBJECT $SESSION"
cd $HOMEDIR/dicoms/$SUBJECT/
tar --remove-files -czf $HOMEDIR/dicoms/$SUBJECT/${SUFFIX}.tar.gz ${SUFFIX} 
rm -r $HOMEDIR/data_BIDS/sourcedata/sub-${SUBJECT}/ses-${SESSION}
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

