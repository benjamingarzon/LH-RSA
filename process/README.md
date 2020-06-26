
# Processing the data:

## Converting the data
process/conversions.sh

## Prepare masks
process/extract_labels.sh

## Functional
process/run_pipeline.sh
  -> process/pipeline.sh

###MPRAGE calculation  
  -> process/mprageconvert/create_mp2rage_command.m

### T1 calculation
  -> process/T1values/MP2RAGE_B1corr/estimateT1.m
### PCM and MVPA calculation
  -> process/setup_GLM.sh
  -> process/prepare_queryengine.py
  -> process/compute_effects.py
  -> process/surface_searchlight.py
      -> utils.py

## Structural freesurfer
process/run_structural_pipeline.sh 
  -> process/run_fs.sh
process/prepare_fs_analysis.sh # prepare thickness and T1 maps for freesurfer analysis

## Structural vbm


# Analyzing the data:
# Structural
analysis/structural_analysis.R

# Functional 
analysis/high_level_analysis.R

