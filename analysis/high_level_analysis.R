rm(list = ls())
library(pracma)
library(ggplot2)
library(dplyr)
library(car)

source('~/Software/ImageLMMR/ImageLMMR.R')
setwd("~/Software/LeftHand/analysis")
source('./high_level_analysis_funcs.R')
source('./mytests.R')
source('./myfmritests.R')
source('./load_covariates.R')

nifti_convert = '/data/lv0/Software/workbench/bin_rh_linux64/wb_command -metric-convert -to-nifti'


# pes : Fixation (1), Stretch (3), TrainedCorrect (5), UntrainedCorrect (7), TrainedIncorrect (9), UntrainedIncorrect (11)
# cope : TrainedCorrect (1), UntrainedCorrect (2), TrainedIncorrect (3), UntrainedIncorrect (4)

###############
# change back : exclude, select first session
###############

WD = '/data/lv0/MotorSkill/fmriprep/analysis'

# define surfaces
rh.gii = '/data/lv0/MotorSkill/fmriprep/freesurfer/fsaverage6/surf/rh.white.surf.gii' 
lh.gii = '/data/lv0/MotorSkill/fmriprep/freesurfer/fsaverage6/surf/lh.white.surf.gii'
mask.rh = '/data/lv0/MotorSkill/labels/fsaverage6/rh.motor.rois.nii.gz'
mask.lh = '/data/lv0/MotorSkill/labels/fsaverage6/lh.motor.rois.nii.gz'
mask_whole = 'mask_whole.nii.gz'   
mask_roi = 'mask.nii.gz'

########################
# Set up analyses
########################
collect_data = F
NPROCS = 8

analysis_type = 'surfL'  #volume, surfR/L 

analysis_name = 'Trained_Untrained'
conditions = c(1, 2)
names(conditions) = c('TrainedCorrect', 'UntrainedCorrect')
peprefix = "/cope"

#analysis_name = 'UntrainedCorrect_UntrainedIncorrect'
#conditions = c(2, 4)
#names(conditions) = c('UntrainedCorrect', 'UntrainedIncorrect')

#analysis_name = 'TrainedCorrect_TrainedIncorrect'
#conditions = c(1, 3)
#names(conditions) = c('TrainedCorrect', 'TrainedIncorrect')


#analysis_name = 'Fixation_Stretch'
#conditions = c(1, 3)
#names(conditions) = c('Fixation', 'Stretch')
#peprefix = "/pe"

########################
# Gather data
########################

# list files, adapt depending on type of analysis
contrasts = NULL
image.list = NULL
condition.list = NULL
for (contrastnum in conditions){
contrast = switch(which(analysis_type == c('volume', 'surfR', 'surfL')),
                  paste0(peprefix, contrastnum, ".nii.gz"),
                  paste0(peprefix, contrastnum, ".func.gii"),
                  paste0(peprefix, contrastnum, ".func.gii"))
contrasts = c(contrasts, contrast)
images.found = system(sprintf("find %s/* | grep %s%s", WD, analysis_type, contrast), intern = TRUE)
image.list = c(image.list, 
               images.found)
condition.list = c(condition.list, 
                   rep(names(conditions)[conditions == contrastnum], length(images.found))
                   )
}

if (analysis_type == 'volume') mysurf = '' else mysurf = ifelse(grepl("R", analysis_type, fixed = T), lh.gii, rh.gii)

if (analysis_type %in% c('volume', 'surfR', 'surfL')) {
  sel = grep(analysis_type, image.list)
  image.list = image.list[sel]
  condition.list = condition.list[sel]}
    
print(image.list)

DATADIR = file.path(WD, "higherlevel", analysis_name)
dir.create(DATADIR)
dir.create(file.path(DATADIR, analysis_type))
dir.create(file.path(DATADIR, analysis_type, "tests"))

setwd(DATADIR)
if (analysis_type == 'volume') {
  setwd(WD)
  
  print("Checking dimensions")
  image.list.clean = mask.list = indices = NULL
  j = 1
  for (imgname in image.list){
    header = check_nifti_header(imgname)
    imgdims = header@dim_[2:4]
    
    if (imgdims[1] != 108 | imgdims[2] != 128 | imgdims[3] != 89) print(paste("Removing", imgname))
    else { 
      image.list.clean = c(image.list.clean, imgname)
      indices = c(indices, j)
      j = j + 1
    }
  }
  condition.list = condition.list[indices]
}

if (collect_data) {
  if (analysis_type == 'volume') {

  masks.found = system("find ../fmriprep/sub-lue*/ses*/func/sub-lue*_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz", intern = TRUE)
  for (imgname in masks.found){
    header = check_nifti_header(imgname)
    imgdims = header@dim_[2:4]
    if (imgdims[1] != 108 | imgdims[2] != 128 | imgdims[3] != 89) print(paste("Removing", imgname))
    else mask.list = c(mask.list, imgname)
  }
   image.list = image.list.clean
   image_list_file = file.path(analysis_type, "image_list.txt")
   write.table(image.list, file = file.path(DATADIR, image_list_file), row.names = F, quote = F, col.names = F)
   
   print("Converting and merging mask files")

   system(sprintf("fslmerge -t all_masks %s", mask.list))
   system("fslmaths all_masks -Tmean -mul `fslnvols all_masks` higherlevel/mask_sum") #`fslnvols      higherlevel/mask_mean all_masks`
   system("fslmaths all_masks -Tmin higherlevel/mask_min")
   system("mri_vol2vol --mov /data/lv0/MotorSkill/cat12/mask.nii.gz --targ higherlevel/mask_min.nii.gz --out higherlevel/GMprob.nii.gz --regheader")
   system(sprintf("fslmaths higherlevel/GMprob -thr 0.1 -bin -mul higherlevel/mask_min higherlevel/mask_min; mv higherlevel/mask_min.nii.gz higherlevel/%s/volume/", analysis_name))
   system(sprintf("fslmaths higherlevel/GMprob -thr 0.1 -bin -mul higherlevel/mask_sum higherlevel/mask_sum; mv higherlevel/mask_sum.nii.gz higherlevel/%s/volume/", analysis_name))
   system("rm all_masks.nii.gz higherlevel/GMprob.nii.gz")

   print("Converting and merging signal files")
   setwd(file.path(DATADIR, "volume"))
   system("rm mask.nii.gz; ln -s mask_min.nii.gz mask.nii.gz")
   system("rm mask_whole.nii.gz; ln -s mask_min.nii.gz mask_whole.nii.gz")
   system("fslmerge -t images `cat image_list.txt`")
   system("flirt -in mask -ref mask -out mask_4mm -interp nearestneighbour -applyisoxfm 4")
   system("flirt -in images -ref mask_4mm -out images_4mm -interp nearestneighbour -applyisoxfm 4")
  } else {
    setwd(DATADIR)
    image_list_file = file.path(analysis_type, "image_list.txt")
    write.table(image.list, file = file.path(DATADIR, image_list_file), row.names = F, quote = F, col.names = F)
    print("Converting and merging files")
    for (filename in image.list) {
      #print(filename)
      system(paste(nifti_convert, 
      filename, 
      paste0(filename, ".nii.gz")
      ))
      
      if (! analysis_type %in% c('volume', 'surfR', 'surfL')) system(sprintf("fslroi %s %s 0 1", filename, filename))
      
    }
    
    
    system(sprintf("IM=`cat %s | sed 's/.gii/.gii.nii.gz/g'`; fslmerge -t %s/images $IM; rm $IM", 
                   image_list_file, analysis_type, analysis_type, analysis_type))
    
    system(sprintf("fslmaths %s/images -Tstd -bin %s/mask_whole", 
                   analysis_type, analysis_type))
    
    system(sprintf('fslmaths %s %s/mask.nii.gz -odt float', 
                   ifelse(analysis_type =='surfR', mask.rh, mask.lh), 
                   analysis_type))
    
    # Convert to gifti for checks
    system(paste(gifti_convert, 
                   sprintf('%s/images.nii.gz', analysis_type), 
                   ifelse(analysis_type =='surfR', rh.gii, lh.gii), 
                   sprintf('%s/images.func.gii', analysis_type)
      ))
    system(sprintf('mri_convert %s/images.func.gii %s/images.func.gii', 
                   analysis_type, analysis_type))  
  }
}



# check results: freeview -f /usr/local/freesurfer/7.1.1-1/subjects/fsaverage6/surf/lh.inflated:overlay=INTERCEPT_coef.func.gii
results.comparison_prereg = doit(file.path(DATADIR, analysis_type),
                                 image.list,
                                 modelcomparisonrun_prereg,
                                 'tests/comparison_prereg',
                                 MASK = mask_roi,
                                 IMAGES_NAME = 'image_list.txt',
                                 IMAGING_NAME = 'images.nii.gz',
                                 conditions = condition.list, 
                                 motion = motion,
                                 to_gifti = mysurf)

results.cubic_prereg = doit(file.path(DATADIR, analysis_type),
                            image.list,
                            testcubicrun_prereg,
                            'tests/cubic_prereg',
                            MASK = mask_roi,
                            IMAGES_NAME = 'image_list.txt',
                            IMAGING_NAME = 'images.nii.gz',
                            conditions = condition.list,
                            motion = motion,
                            to_gifti = mysurf)

results.asymptotic_prereg = doit(file.path(DATADIR, analysis_type),
                             image.list,
                            testasymptoticrun_prereg,
                            'tests/asymptotic_prereg',
                            MASK = mask_roi,
                            IMAGES_NAME = 'image_list.txt',
                            IMAGING_NAME = 'images.nii.gz',
                            conditions = condition.list,
                            motion = motion,
                            to_gifti = mysurf)

results.quadratic_prereg = doit(file.path(DATADIR, analysis_type),
                                 image.list,
                                 testquadraticrun_prereg,
                                 'tests/quadratic_prereg',
                                 MASK = mask_roi,
                                 IMAGES_NAME = 'image_list.txt',
                                 IMAGING_NAME = 'images.nii.gz',
                                 conditions = condition.list,
                                 motion = motion,
                                 to_gifti = mysurf)



results.cubic_whole_prereg = doit(file.path(DATADIR, analysis_type),
                                  image.list,
                                  testcubicrun_prereg,
                                  'tests/cubic_whole_prereg',
                                  MASK = mask_whole,
                                  IMAGES_NAME = 'image_list.txt',
                                  IMAGING_NAME = 'images.nii.gz',
                                  conditions = condition.list,
                                  motion = motion,
                                  to_gifti = mysurf)

results.comparison_whole_prereg = doit(file.path(DATADIR, analysis_type),
                                       image.list,
                                       modelcomparisonrun_prereg,
                                       'tests/comparison_whole_prereg_cond',
                                       MASK = mask_whole,
                                       IMAGES_NAME = 'image_list.txt',
                                       IMAGING_NAME = 'images.nii.gz',
                                       conditions = condition.list, 
                                       motion = motion,
                                       to_gifti = mysurf)


stophere





results.asymptotic_prereg_groupxtraining = doit(file.path(DATADIR, analysis_type),
                                                image.list,
                                                testasymptoticrun_prereg_groupxtraining,
                                                'tests/asymptotic_prereg_groupxtraining',
                                                MASK = mask_roi,
                                                IMAGES_NAME = 'image_list.txt',
                                                IMAGING_NAME = 'images.nii.gz',
                                                conditions = condition.list,
                                                motion = motion,
                                                to_gifti = mysurf)

results.asymptotic_prereg_omni = doit(file.path(DATADIR, analysis_type),
                                                image.list,
                                                testasymptoticrun_prereg_omni,
                                                'tests/asymptotic_prereg_omni',
                                                MASK = mask_roi,
                                                IMAGES_NAME = 'image_list.txt',
                                                IMAGING_NAME = 'images.nii.gz',
                                                conditions = condition.list,
                                                motion = motion,
                                                to_gifti = mysurf)

stophere
results.testgeneralization_prereg = doit(file.path(DATADIR, analysis_type),
                                       image.list,
                                       testgeneralization_prereg,
                                       'tests/generalization_prereg',
                                       MASK = mask_roi,
                                       IMAGES_NAME = 'image_list.txt',
                                       IMAGING_NAME = 'images.nii.gz',
                                       conditions = condition.list,
                                       motion = motion,
                                       to_gifti = mysurf)


stophere
results.testfirstsessionaverage = doit(file.path(DATADIR, analysis_type),
                                      image.list,
                                      testfirstsessionaverage_prereg,
                                      'tests/firstsessionaverage_prereg',
                                      MASK = mask_whole,
                                      IMAGES_NAME = 'image_list.txt',
                                      IMAGING_NAME = 'images.nii.gz',
                                      conditions = condition.list,
                                      motion = motion,
                                      to_gifti = mysurf)


results.asymptotic_whole_prereg = doit(file.path(DATADIR, analysis_type),
                                       image.list,
                                       testasymptoticrun_prereg,
                                       'tests/asymptotic_whole_prereg',
                                       MASK = mask_whole,
                                       IMAGES_NAME = 'image_list.txt',
                                       IMAGING_NAME = 'images.nii.gz',
                                       conditions = condition.list,
                                       motion = motion,
                                       to_gifti = mysurf)


results.quadratic_whole_prereg = doit(file.path(DATADIR, analysis_type),
                                      image.list,
                                      testquadraticrun_prereg,
                                      'tests/quadratic_whole_prereg',
                                      MASK = mask_whole,
                                      IMAGES_NAME = 'image_list.txt',
                                      IMAGING_NAME = 'images.nii.gz',
                                      conditions = condition.list,
                                      motion = motion,
                                      to_gifti = mysurf)
stophere
# tests for activation maps
results.asymptotic.prereg = doit(file.path(DATADIR, analysis_type),
                               image.list,
                               testasymptoticrun_prereg,
                               'tests/asymptotic_prereg',
                               MASK = mask_roi,
                               IMAGES_NAME = 'image_list.txt',
                               IMAGING_NAME = 'images.nii.gz',
                               conditions = condition.list,
                               motion = motion,
                               to_gifti = mysurf)


results.comparison.prereg = doit(file.path(DATADIR, analysis_type),
                                 image.list,
                                 modelcomparisonrun,
                                 'tests/comparison_prereg',
                                 MASK = mask_roi,
                                 IMAGES_NAME = 'image_list.txt',
                                 IMAGING_NAME = 'images.nii.gz',
                                 conditions = condition.list, 
                                 motion = motion,
                                 to_gifti = mysurf)



results.quadratic_prereg_groupxtraining = doit(file.path(DATADIR, analysis_type),
                                image.list,
                                testquadraticrun_prereg_groupxtraining,
                                'tests/quadratic_prereg_groupxtraining',
                                MASK = mask_roi,
                                IMAGES_NAME = 'image_list.txt',
                                IMAGING_NAME = 'images.nii.gz',
                                conditions = condition.list,
                                motion = motion,
                                to_gifti = mysurf)

results.asymptotic_prereg_groupxtraining = doit(file.path(DATADIR, analysis_type),
                                               image.list,
                                               testasymptoticrun_prereg_groupxtraining,
                                               'tests/asymptotic_prereg_groupxtraining',
                                               MASK = mask_roi,
                                               IMAGES_NAME = 'image_list.txt',
                                               IMAGING_NAME = 'images.nii.gz',
                                               conditions = condition.list,
                                               motion = motion,
                                               to_gifti = mysurf)

results.quadratic_prereg_groupxconditionxtraining = doit(file.path(DATADIR, analysis_type),
                                               image.list,
                                               testquadraticrun_prereg_groupxconditionxtraining,
                                               'tests/quadratic_prereg_groupxconditionxtraining',
                                               MASK = mask_roi,
                                               IMAGES_NAME = 'image_list.txt',
                                               IMAGING_NAME = 'images.nii.gz',
                                               conditions = condition.list,
                                               motion = motion,
                                               to_gifti = mysurf)

results.asymptotic_prereg_groupxconditionxtraining = doit(file.path(DATADIR, analysis_type),
                                                         image.list,
                                                         testasymptoticrun_prereg_groupxconditionxtraining,
                                                         'tests/asymptotic_prereg_groupxconditionxtraining',
                                                         MASK = mask_roi,
                                                         IMAGES_NAME = 'image_list.txt',
                                                         IMAGING_NAME = 'images.nii.gz',
                                                         conditions = condition.list,
                                                         motion = motion,
                                                         to_gifti = mysurf)

results.quadratic.prereg = doit(file.path(DATADIR, analysis_type),
                                image.list,
                                testquadraticrun_prereg,
                                'tests/quadratic_prereg',
                                MASK = mask_roi,
                                IMAGES_NAME = 'image_list.txt',
                                IMAGING_NAME = 'images.nii.gz',
                                conditions = condition.list,
                                motion = motion,
                                to_gifti = mysurf)

results.asymptotic.prereg = doit(file.path(DATADIR, analysis_type),
                                image.list,
                                testasymptoticrun_prereg,
                                'tests/asymptotic_prereg',
                                MASK = mask_roi,
                                IMAGES_NAME = 'image_list.txt',
                                IMAGING_NAME = 'images.nii.gz',
                                conditions = condition.list,
                                motion = motion,
                                to_gifti = mysurf)



results.linear.prereg = doit(file.path(DATADIR, analysis_type),
                             image.list,
                             testlinearrun_prereg,
                             'tests/linear_prereg',
                             MASK = mask_roi,
                             IMAGES_NAME = 'image_list.txt',
                             IMAGING_NAME = 'images.nii.gz',
                             conditions = condition.list,
                             motion = motion,
                             to_gifti = mysurf)

  if (F) {
    
results.linear = doit(file.path(DATADIR, analysis_type),
                           image.list,
                           testlinearrun,
                           'tests/linear',
                           MASK = mask_roi,
                           IMAGES_NAME = 'image_list.txt',
                           IMAGING_NAME = 'images.nii.gz',
                           conditions = condition.list,
                           motion = motion,
                           to_gifti = mysurf)


results.comparison_whole_prereg = doit(file.path(DATADIR, analysis_type),
                          image.list,
                          modelcomparisonrun,
                          'tests/comparison_whole_prereg',
                          MASK = mask_whole,
                          IMAGES_NAME = 'image_list.txt',
                          IMAGING_NAME = 'images.nii.gz',
                          conditions = condition.list, 
                          motion = motion,
                          to_gifti = mysurf)

results.comparison.prereg = doit(file.path(DATADIR, analysis_type),
                                 image.list,
                                 modelcomparisonrun,
                                 'tests/comparison_prereg',
                                 MASK = mask_roi,
                                 IMAGES_NAME = 'image_list.txt',
                                 IMAGING_NAME = 'images.nii.gz',
                                 conditions = condition.list, 
                                 motion = motion,
                                 to_gifti = mysurf)

results.linear_whole_prereg = doit(file.path(DATADIR, analysis_type),
                                   image.list,
                                   testlinearrun_prereg,
                                   'tests/linear_whole_prereg',
                                   MASK = mask_whole,
                                   IMAGES_NAME = 'image_list.txt',
                                   IMAGING_NAME = 'images.nii.gz',
                                   conditions = condition.list,
                                   motion = motion,
                                   to_gifti = mysurf)


results.linear.prereg_half = doit(file.path(DATADIR, analysis_type),
                                  image.list,
                                  testlinearrun_prereg_half,
                                  'tests/linear_prereg_half',
                                  MASK = mask_roi,
                                  IMAGES_NAME = 'image_list.txt',
                                  IMAGING_NAME = 'images.nii.gz',
                                  conditions = condition.list,
                                  motion = motion,
                                  to_gifti = mysurf)

}

# sync folders
#mydir=Trained_Untrained
#for name in comparison_prereg comparison_whole_prereg linear_prereg linear linear_whole_prereg linear_prereg_half quadratic_prereg quadratic_whole_prereg quadratic_prereg_groupxtraining quadratic_prereg_groupxconditionxtraining;# 
#for name in asymptotic_prereg_groupxtraining asymptotic_prereg_groupxconditionxtraining;# 
#do 
#ln -s /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surfL/tests/$name /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surf/tests/${name}.lh
#ln -s /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surfR/tests/$name /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surf/tests/${name}.rh
#done


# mydir=Trained_Untrained
# for name in comparison_prereg_ext comparison_whole_prereg_ext cubic_prereg  cubic_whole_prereg; do 
# ln -s /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surfL/tests/$name /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surf/tests/${name}.lh
# ln -s /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surfR/tests/$name /data/lv0/MotorSkill/fmriprep/analysis/higherlevel/$mydir/surf/tests/${name}.rh
# done
