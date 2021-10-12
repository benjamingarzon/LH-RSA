#!/bin/sh
#suffix='tess' # whole brain using tessellation
#labels_file='/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/mask_tessellation162_parcels.txt'
#labels_dir='/data/lv0/MotorSkill/labels/fsaverage/tessellation162'
#base_suffix=tess-cross

labels_file='/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/motor_roi_parcels.txt'
labels_dir='/data/lv0/MotorSkill/labels/fsaverage'
base_suffix=mask-cross

WD='/data/lv0/MotorSkill/'
num_cores=10
n_sample=100

# cross-validated
python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=run \
    --suffix=${base_suffix}-runprew \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=effects.nii.gz \
    --n_sample=$n_sample 
    
python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=run \
    --suffix=${base_suffix}-runprew-perm \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=effects.nii.gz \
    --n_sample=$n_sample \
    --permutate  &

#    --just_gather \

exit 1

python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --suffix=${base_suffix}-noprew \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=effects.nii.gz \
    --n_sample=$n_sample  &
    
# cross-validated & permutated
python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=run \
    --suffix=${base_suffix}-runprew-perm \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=effects.nii.gz \
    --n_sample=$n_sample \
    --permutate  &

python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --suffix=${base_suffix}-noprew-perm \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=effects.nii.gz \
    --n_sample=$n_sample \
    --permutate 

exit 1
# cross-validated w derivatives & permutated
python surface_roi_analysis.py --WD=$WD \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=session \
    --suffix=${base_suffix}-derivatives \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=derivatives.nii.gz \
    --permutate &

exit 1
    
exit 1
# cross-validated w derivatives
python surface_roi_analysis.py --WD=$WD \
    --overwrite_extract \
    --overwrite_scores \
    --output_data \
    --do_prewhitening=session \
    --suffix=${base_suffix}-derivatives \
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
    --suffix=${base_suffix}-derivatives \
    --num_cores=$num_cores \
    --labels_file=$labels_file \
    --labels_dir=$labels_dir \
    --effects_name=derivatives.nii.gz \
    --permutate &
