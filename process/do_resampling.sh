#/bin/sh
# do resampling

SUBJECT=$1
SESS=$2
RUN1=1
RUN2=2
MYDIR=/home/xgarzb@GU.GU.SE/Data/LeftHand/Lund1/fmriprep/analysis/$SUBJECT/ses-$SESS/

cd $MYDIR

if [ ]; then
mv run${RUN1}/effects_LSA_${RUN1}.nii.gz run${RUN1}/effects_LSA_${RUN1}_orig.nii.gz
mri_vol2vol --mov  run${RUN1}/effects_LSA_${RUN1}_orig.nii.gz  --targ run${RUN2}/effects_LSA_${RUN2}.nii.gz --o run${RUN1}/effects_LSA_${RUN1}.nii.gz --regheader

fslmerge -t effects `ls -v run*/effects_LSA_?.nii.gz` 
cd ../..
fi
for i in 1 2 3 4; do
  
  mv run${RUN1}/volume/cope${i}.nii.gz run${RUN1}/volume/cope${i}_orig.nii.gz
  mri_vol2vol --mov run${RUN1}/volume/cope${i}_orig.nii.gz --targ run${RUN2}/volume/cope${i}.nii.gz --o run${RUN1}/volume/cope${i}.nii.gz --regheader
done

#for image in ; do
#mv run${RUN1}/effects_LSA_${RUN1}.nii.gz run${RUN1}/effects_LSA_${RUN1}_orig.nii.gz
#done

#for i in `cat image_list.txt`; do val=`fslval $i dim1`; if [ ! $val == 108 ]; then echo $i $val; fi; done
#for i in `find sub-lue*/ses*/run*/volume/cope1.nii.gz`; do val=`fslval $i dim1`; if [ ! $val == 108 ]; then echo $i $val; fi; done
# for i in `find */sub-lue*/ses*/func/sub-lue*_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz`; do val=`fslval $i dim1`; if [ ! $val == 108 ]; then echo $i $val; fi; done
#./do_resampling.sh sub-lue2102 3 #
#./do_resampling.sh sub-lue2105 5 #
#./do_resampling.sh sub-lue3104 2 #
#./do_resampling.sh sub-lue3205 2 #
#./do_resampling.sh sub-lue3205 3 #
#/do_resampling.sh sub-lue2106 6
#./do_resampling.sh sub-lue3106 1
#./do_resampling.sh sub-lue4205 1
#./do_resampling.sh sub-lue5104 1


