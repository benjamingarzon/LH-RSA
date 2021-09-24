#!/bin/sh
#suffix='tess' # whole brain using tessellation
#labels_file='/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/tessellation162_parcels.txt'
#labels_dir='/data/lv0/MotorSkill/labels/fsaverage/tessellation162'

labels_file='/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/motor_roi_parcels.txt'
labels_dir='/data/lv0/MotorSkill/labels/fsaverage'
WD='/data/lv0/MotorSkill/'
num_cores=10
# cross-validated & permutated
python surface_roi_analysis.py --WD=$WD \
    --output_data \
    --do_prewhitening=session \
    --suffix=mask-cross-perm \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=effects.nii.gz \
    --permutate &
exit 1
# cross-validated w derivatives & permutated
python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=session \
    --suffix=mask-cross-derivatives \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=derivatives.nii.gz \
    --permutate &

exit 1
# cross-validated
python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=session \
    --suffix=mask-cross \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=effects.nii.gz &
    
exit 1
# cross-validated w derivatives
python surface_roi_analysis.py --WD=$WD \
    --overwrite_extract \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=session \
    --suffix=mask-cross-derivatives \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=derivatives.nii.gz &

# cross-validated w derivatives & permutated
python surface_roi_analysis.py --WD=$WD \
    --overwrite_extract \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=session \
    --suffix=mask-cross-derivatives \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=derivatives.nii.gz \
    --permutate &
