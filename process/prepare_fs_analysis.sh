#!/usr/bin/sh
# different metrics: area, t1 values
HOMEDIR=/home/benjamin.garzon/Data/LeftHand/Lund1
STRUCT_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
PROGDIR=~/Software/LeftHand/process/
T1DIR=/home/benjamin.garzon/Software/LeftHand/process/T1values/MP2RAGE_B1corr
SUBJECTS_DIR=/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer
SUBJECTS_DIR=/home/share/MotorSkill/freesurfer/

rm /home/share/MotorSkill/freesurfer/fsaverage
ln -s /usr/local/freesurfer/subjects/fsaverage /home/share/MotorSkill/freesurfer/fsaverage
cd $SUBJECTS_DIR
SUFFIX=$3
SUBJECTS=`echo sub-lue${SUFFIX}??? | sed 's/sub-//g'`
echo $SUBJECTS

if [ $1 ]; then TYPE=$1; else TYPE=thickness; fi
if [ $2 ]; then FRAC=$2; else FRAC=0.25; fi

HCPDIR=/home/share/Software/HCP/workbench/bin_rh_linux64/

WD=/home/benjamin.garzon/Data/LeftHand/Lund1/data_BIDS
NSESSIONS=7
SMOOTH=10
OVERWRITE="0"
for SUBJECT in $SUBJECTS; do
    for SESSION in `seq $NSESSIONS`; do
       for HEMI in rh lh; do           
           
	   if [ $TYPE == "T1" ]; then
               METRIC="${TYPE}_${FRAC}"

  	       if [ "$OVERWRITE" == "0" ] && \
		      [ -e "$SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.nii.gz" ]; then
	 		echo "Already done: `ls -lh $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.nii.gz`" 
			continue
	       fi

     	       if [ ! -e "$SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness" ]; then
			echo "Not found: $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.thickness"
			continue
	       fi


	       if [ ! -e $HOMEDIR/vbm/sub-${SUBJECT}/ses-${SESSION}/T1map.nii.gz ]; then 
	           matlab -nosplash -nodisplay -r "addpath('$T1DIR'); estimateT1('$HOMEDIR/vbm/sub-${SUBJECT}/ses-${SESSION}/T1wtotemplate_brain.nii.gz', '$HOMEDIR/vbm/sub-${SUBJECT}/ses-${SESSION}/T1map.nii.gz'); exit"
	       fi

	       cd $PROGDIR/T1values/MP2RAGE_B1corr
               mri_vol2surf --mov $HOMEDIR/vbm/sub-${SUBJECT}/ses-${SESSION}/T1map.nii.gz \
               --regheader  sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base --surf white \
               --projfrac $FRAC --surf-fwhm $SMOOTH --interp trilinear --cortex --hemi $HEMI --out_type mgh \
               --o $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.mgh

                mris_convert -c $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.mgh \
                $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white \
                $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.gii
           else
             	METRIC=$TYPE
    	        if [ "$OVERWRITE" == "0" ] && \
		      [ -e "$SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.nii.gz" ]; then
	 		echo "Already done: `ls -lh $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.nii.gz`" 
			continue
	        fi

		if [ ! -e "$SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}" ]; then
			echo "Not found: $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}"
			continue
		fi

	        mris_convert -c $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC} \
	        $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white \
	        $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.gii
	   fi	

    
           mris_convert $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii
    
           mris_convert $SUBJECTS_DIR/fsaverage/surf/${HEMI}.sphere.reg \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii

           mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/${HEMI}.white \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white.fsaverage.gii
           
           $HCPDIR/wb_command -metric-resample \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii \
           BARYCENTRIC \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.func.gii

           $HCPDIR/wb_command -metric-smoothing \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white.fsaverage.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.func.gii \
           $SMOOTH \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.$SMOOTH.func.gii
                          
           $HCPDIR/wb_command -metric-convert -to-nifti \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.$SMOOTH.func.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.nii.gz 
           
           rm $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.${METRIC}.*gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.sphere.reg.fsaverage.gii \
           $SUBJECTS_DIR/sub-${SUBJECT}.${SESSION}.long.sub-${SUBJECT}.base/surf/${HEMI}.white.fsaverage.gii

       done # hemi

    done
    
    # gather volumes 
    fslmerge -t $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/rh.${METRIC}.nii.gz $SUBJECTS_DIR/sub-${SUBJECT}.?.long.sub-${SUBJECT}.base/surf/rh.${METRIC}.nii.gz
    fslmerge -t $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/lh.${METRIC}.nii.gz $SUBJECTS_DIR/sub-${SUBJECT}.?.long.sub-${SUBJECT}.base/surf/lh.${METRIC}.nii.gz
    fslmaths $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/rh.${METRIC}.nii.gz -Tstd $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/rh.${METRIC}.std.nii.gz
    fslmaths $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/lh.${METRIC}.nii.gz -Tstd $SUBJECTS_DIR/sub-${SUBJECT}.base/surf/lh.${METRIC}.std.nii.gz

done

mkdir $SUBJECTS_DIR/results/
# gather all
ls $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/rh.${METRIC}.nii.gz | cut -d'-' -f2 | cut -d '.' --output-delimiter ' ' -f1,2 > $SUBJECTS_DIR/results/rh.${METRIC}.txt
ls $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/lh.${METRIC}.nii.gz | cut -d'-' -f2 | cut -d '.' --output-delimiter ' ' -f1,2 > $SUBJECTS_DIR/results/lh.${METRIC}.txt

fslmerge -t $SUBJECTS_DIR/results/rh.${METRIC}.nii.gz $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/rh.${METRIC}.nii.gz 
fslmerge -t $SUBJECTS_DIR/results/lh.${METRIC}.nii.gz $SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/lh.${METRIC}.nii.gz 

mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/lh.white \
$SUBJECTS_DIR/results/lh.white.gii

mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/rh.white \
$SUBJECTS_DIR/results/rh.white.gii

fslmaths $SUBJECTS_DIR/results/lh.${METRIC}.nii.gz -Tmean $SUBJECTS_DIR/results/lh.${METRIC}_mean.nii.gz
fslmaths $SUBJECTS_DIR/results/rh.${METRIC}.nii.gz -Tmean $SUBJECTS_DIR/results/rh.${METRIC}_mean.nii.gz

$HCPDIR/wb_command -metric-convert -from-nifti \
$SUBJECTS_DIR/results/lh.${METRIC}_mean.nii.gz \
$SUBJECTS_DIR/results/lh.white.gii \
$SUBJECTS_DIR/results/lh.${METRIC}_mean.func.gii

$HCPDIR/wb_command -metric-convert -from-nifti \
$SUBJECTS_DIR/results/rh.${METRIC}_mean.nii.gz \
$SUBJECTS_DIR/results/rh.white.gii \
$SUBJECTS_DIR/results/rh.${METRIC}_mean.func.gii

rm $SUBJECTS_DIR/results/?h.${METRIC}_mean.nii.gz 

# mean of std
if [ ]; then
for GROUP in 1 2; do 
    fslmerge -t $SUBJECTS_DIR/results/rh.${METRIC}.${GROUP}.std.nii.gz $SUBJECTS_DIR/sub-lue1${GROUP}*.base/surf/rh.${METRIC}.std.nii.gz
    fslmerge -t $SUBJECTS_DIR/results/lh.${METRIC}.${GROUP}.std.nii.gz $SUBJECTS_DIR/sub-lue1${GROUP}*.base/surf/lh.${METRIC}.std.nii.gz
    fslmaths $SUBJECTS_DIR/results/rh.${METRIC}.${GROUP}.std.nii.gz -Tmean $SUBJECTS_DIR/results/rh.${METRIC}.${GROUP}.std.nii.gz
    fslmaths $SUBJECTS_DIR/results/lh.${METRIC}.${GROUP}.std.nii.gz -Tmean $SUBJECTS_DIR/results/lh.${METRIC}.${GROUP}.std.nii.gz
    
    $HCPDIR/wb_command -metric-convert -from-nifti \
    $SUBJECTS_DIR/results/lh.${METRIC}.${GROUP}.std.nii.gz \
    $SUBJECTS_DIR/results/lh.white.gii \
    $SUBJECTS_DIR/results/lh.${METRIC}.${GROUP}.std.func.gii
    
    $HCPDIR/wb_command -metric-convert -from-nifti \
    $SUBJECTS_DIR/results/rh.${METRIC}.${GROUP}.std.nii.gz \
    $SUBJECTS_DIR/results/rh.white.gii \
    $SUBJECTS_DIR/results/rh.${METRIC}.${GROUP}.std.func.gii

done 
fi
rm $SUBJECTS_DIR/results/?h.${METRIC}.*std.nii.gz $SUBJECTS_DIR/results/?h.${METRIC}_mean.nii.gz #$SUBJECTS_DIR/sub-*.?.long.sub-*.base/surf/?h.${METRIC}.nii.gz
rm $SUBJECTS_DIR/sub-*.base/surf/?h.${METRIC}.?.std.nii.gz $SUBJECTS_DIR/results/?h.white.gii

cd $SUBJECTS_DIR/results

##mris_calc -o rh.${METRIC}.std.diff.func.gii rh.${METRIC}.1.std.func.gii sub rh.${METRIC}.2.std.func.gii
##mris_calc -o lh.${METRIC}.std.diff.func.gii lh.${METRIC}.1.std.func.gii sub lh.${METRIC}.2.std.func.gii

#freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:overlay=$SUBJECTS_DIR/results/lh.${METRIC}.std.diff.func.gii $SUBJECTS_DIR/fsaverage/surf/rh.inflated:overlay=$SUBJECTS_DIR/results/rh.${METRIC}.std.diff.func.gii
#freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:overlay=$SUBJECTS_DIR/results/lh.${METRIC}_mean.func.gii $SUBJECTS_DIR/fsaverage/surf/rh.inflated:overlay=$SUBJECTS_DIR/results/rh.${METRIC}_mean.func.gii

#freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:overlay=lh.${METRIC}.1.std.func.gii $SUBJECTS_DIR/fsaverage/surf/rh.inflated:overlay=rh.${METRIC}.1.std.func.gii
#freeview -f f$SUBJECTS_DIR/saverage/surf/lh.inflated:overlay=lh.${METRIC}.2.std.func.gii $SUBJECTS_DIR/fsaverage/surf/rh.inflated:overlay=rh.${METRIC}.2.std.func.gii

##fslmaths $SUBJECTS_DIR/results/lh.${METRIC}.nii.gz -Tmean -bin $SUBJECTS_DIR/results/lh.${METRIC}.mask.nii.gz
##fslmaths $SUBJECTS_DIR/results/rh.${METRIC}.nii.gz -Tmean -bin $SUBJECTS_DIR/results/rh.${METRIC}.mask.nii.gz

# create a mask in motor cortex
#fslmaths $SUBJECTS_DIR/results/lh.${METRIC}.nii.gz -Tmean -thr 3 -bin $SUBJECTS_DIR/results/lh.${METRIC}.small.mask.nii.gz
#fslmaths $SUBJECTS_DIR/results/rh.${METRIC}.nii.gz -Tmean -thr 3 -bin $SUBJECTS_DIR/results/rh.${METRIC}.smallmask.nii.gz

mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/lh.pial $SUBJECTS_DIR/results/lh.fsaverage.white.gii
mris_convert --to-scanner $SUBJECTS_DIR/fsaverage/surf/rh.pial $SUBJECTS_DIR/results/rh.fsaverage.white.gii

# if there are several T1 depths merge together in 1 file
if [ $TYPE == "T1" ]; then
fslmerge -t lh.T1.nii.gz lh.T1_*.nii.gz
fslmerge -t rh.T1.nii.gz rh.T1_*.nii.gz
rm lh.T1.txt rh.T1.txt
for d in $DEPTHS; do
echo $d
sed "s/$/\t$d/" lh.T1_${d}.txt >> lh.T1.txt
sed "s/$/\t$d/" rh.T1_${d}.txt >> rh.T1.txt
done
fi

cp 
