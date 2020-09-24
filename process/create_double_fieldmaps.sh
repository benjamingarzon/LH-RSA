#/bin/sh
# in 3 datasets there were double fieldmaps
# xx lue5104       5  US_20200223_15.26.09     601  1801    4,5
# xx lue4207       4  US_20191208_11.02.31     601  1701  3,4,5
# xx lue5202       7  US_20200307_15.48.28     701  1401  3,4,5

BIDSDIR=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
DATADIR=/home/share/MotorSkill
SUBJECT=lue4207
SESSION=4
WAVE=4
USDIR=US_20191208_11.02.31
B0DIR="1701_B0map"
RUNS1="1 2"
RUNS2="3 4 5"
FMAPDIR=$BIDSDIR/sub-$SUBJECT/ses-$SESSION/fmap

ls $FMAPDIR

# backup just in case
cp -r $FMAPDIR ${FMAPDIR}_orig

# uncompress
cd $DATADIR/dicoms_wave$WAVE/$SUBJECT/ 
tar -xvzf $DATADIR/dicoms_wave$WAVE/$SUBJECT/${USDIR}.tar.gz

# convert and rename
cd $DATADIR/dicoms_wave$WAVE/$SUBJECT/${USDIR}/Dicom/
cd ${B0DIR}*
dcm2niix *
gzip *.nii

#phase
mv *e1.nii.gz sub-${SUBJECT}_ses-${SESSION}_epi2.nii.gz
mv *e1.json sub-${SUBJECT}_ses-${SESSION}_epi2.json
mv *e1a.nii.gz sub-${SUBJECT}_ses-${SESSION}_epi1.nii.gz
mv *e1a.json sub-${SUBJECT}_ses-${SESSION}_epi1.json

#assign runs
for RUN in $RUNS2; do
cp sub-${SUBJECT}_ses-${SESSION}_epi2.nii.gz sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi2.nii.gz
cp sub-${SUBJECT}_ses-${SESSION}_epi2.json sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi2.json
cp sub-${SUBJECT}_ses-${SESSION}_epi1.nii.gz sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi1.nii.gz
cp sub-${SUBJECT}_ses-${SESSION}_epi1.json sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi1.json
done
rm sub-${SUBJECT}_ses-${SESSION}_epi?.*

mv *.json ${FMAPDIR}
mv *.nii.gz ${FMAPDIR}

# assign runs

cd $FMAPDIR
for RUN in $RUNS1; do
cp sub-${SUBJECT}_ses-${SESSION}_epi2.nii.gz sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi2.nii.gz
cp sub-${SUBJECT}_ses-${SESSION}_epi2.json sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi2.json
cp sub-${SUBJECT}_ses-${SESSION}_epi1.nii.gz sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi1.nii.gz
cp sub-${SUBJECT}_ses-${SESSION}_epi1.json sub-${SUBJECT}_ses-${SESSION}_run-0${RUN}_epi1.json
done

rm sub-${SUBJECT}_ses-${SESSION}_epi?.* *fieldmap* magnitude*
rm -r $DATADIR/dicoms_wave$WAVE/$SUBJECT/${USDIR}