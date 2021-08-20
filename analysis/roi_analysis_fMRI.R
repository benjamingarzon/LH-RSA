

# load model
rm(list = ls())
cat("\014")

# plot data

# run only for Classic within mask
# plot t1 from where there is an effect
# write how many observations
# add annotation :annot=/data/lv0/MotorSkill/labels/subject/fsaverage6/label/rh.motor.rois.annot:annot_outline=1


setwd("~/Software/LeftHand/analysis")
source("./plot_funcs.R")

NCOL = 3

#FIGS_DIR = '~/Data/LeftHand/Lund1/figs/Trained_Untrained/'
FIGS_DIR = '~/Data/LeftHand/Lund1/figs/TrainedCorrect_TrainedIncorrect/'
WD = '/data/lv0/MotorSkill/fmriprep/analysis/higherlevel/TrainedCorrect_TrainedIncorrect/'
dir.create(FIGS_DIR)
#unlink(file.path(FIGS_DIR, 'figs','*.png'))

THR = 0.99

#################################################################################
radius = 3

surface.rh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/rh.pial.shape.gii'
surface.lh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/lh.pial.shape.gii'
SURFACE_FILES = c(surface.lh, surface.rh)
names(SURFACE_FILES) = c("lh", "rh")

inflated.rh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/rh.inflated.shape.gii:curvature=/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/rh.curv'
inflated.lh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/lh.inflated.shape.gii:curvature=/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/lh.curv'

DISTANCE = 10

#WD = '/data/lv0/MotorSkill/fmriprep/analysis/higherlevel/TrainedCorrect_TrainedIncorrect/'
MASK_NAME = 'mask.nii.gz'

if (F) {
DATADIR = file.path(WD, 'volume/')

TESTDIR = 'tests/linear'
TESTNAME = 'GROUPExp_x_CONDITIONUnt_x_TRAINING_p' #'TRAINING_p_fdr'
myplots.vol.training = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                       plot_function = plot_activation_data)
TESTDIR = 'tests/linear'
TESTNAME = 'CONDITIONUnt_p' #'TRAINING_p_fdr'
myplots.vol.training = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                       plot_function = plot_activation_data)



DATADIR = file.path(WD, 'volume/')
TESTDIR = 'tests/linear'
TESTNAME = 'GROUPExp_x_CONDITIONUnt_x_TRAINING_p' #'TRAINING_p_fdr'
myplots.vol.training = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                    plot_function = plot_activation_data)

}

# surf
radius = 10
annot=list('lh'='/data/lv0/MotorSkill/labels/subject/fsaverage6/label/lh.motor.rois.annot',
           'rh'='/data/lv0/MotorSkill/labels/subject/fsaverage6/label/rh.motor.rois.annot')
aparc=list('lh'='/data/lv0/MotorSkill/labels/subject/fsaverage6/label/lh.aparc.annot',
           'rh'='/data/lv0/MotorSkill/labels/subject/fsaverage6/label/rh.aparc.annot')

THR = 0.975

DATADIR = file.path(WD, 'surf/')

MASK_NAME = 'mask_whole.nii.gz'
TESTDIR = 'tests/quadratic_whole_prereg'
TESTNAME = 'INTERCEPT_p_fdr'
myplots.surf.intercept = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                          plot_function = plot_activation_data, annot = aparc)

stophere
MASK_NAME = 'mask_prereg.nii.gz'
TESTDIR = 'tests/quadratic_prereg'  
TESTNAME = 'Omni_p_fdr'
myplots.surf.omni = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                     plot_function = plot_activation_data, annot = annot)

MASK_NAME = 'mask_prereg.nii.gz'
TESTDIR = 'tests/quadratic_prereg_groupxconditionxtraining'
TESTNAME = 'Omni_p_fdr'
myplots.surf.omni = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                          plot_function = plot_activation_data, annot = annot)


MASK_NAME = 'mask_prereg.nii.gz'
TESTDIR = 'tests/quadratic_prereg_groupxtraining'
TESTNAME = 'Omni_p_fdr'
myplots.surf.omni = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                     plot_function = plot_activation_data, annot = annot)

#MASK_NAME = 'mask_prereg.nii.gz'
#TESTDIR = 'tests/linear_prereg'
#TESTNAME = 'INTERCEPT_p_fdr'
#myplots.surf.intercept = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
#                                          plot_function = plot_activation_data, annot = annot)

MASK_NAME = 'mask_whole.nii.gz'
TESTDIR = 'tests/quadratic_whole_prereg'
TESTNAME = 'Omni_p_fdr'
myplots.surf.omni = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                          plot_function = plot_activation_data, annot = aparc)

MASK_NAME = 'mask_whole.nii.gz'
TESTDIR = 'tests/quadratic_whole_prereg'
TESTNAME = 'INTERCEPT_p_fdr'
myplots.surf.intercept = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                          plot_function = plot_activation_data, annot = aparc)

MASK_NAME = 'mask_whole.nii.gz'
TESTDIR = 'tests/quadratic_whole_prereg'
TESTNAME = 'TRAINING_p_fdr'
myplots.surf.training = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                          plot_function = plot_activation_data, annot = aparc)

stophere

MASK_NAME = 'mask_whole.nii.gz'
TESTDIR = 'tests/linear_whole_prereg'
TESTNAME = 'TRAINING_p'
myplots.surf.training = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                         plot_function = plot_activation_data, annot = aparc)



TESTDIR = 'tests/linear_prereg'  
TESTNAME = 'GROUPExp_x_CONDITIONUnt_p_fdr'
THR = 0.975
myplots.surf.group_x_condition = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                         plot_function = plot_activation_data, annot = annot)



stophere
NCOL = 2
TESTDIR = 'tests/linear'  # no clusters
TESTNAME = 'GROUPExp_x_CONDITIONUnt_x_TRAINING_p'
THR = 0.975
myplots.surf.group_x_training = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                         plot_function = plot_activation_data, annot = annot)

stophere
TESTDIR = 'tests/linear_prereg'
TESTNAME = 'CONDITIONUnt_p_fdr'
THR = 0.95
#myplots.surf.condition = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
#                                          plot_function = plot_activation_data, annot = annot)

TESTDIR = 'tests/linear_prereg'
TESTNAME = 'TRAINING_p_fdr'
THR = 0.95
#myplots.surf.training = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
#                                         plot_function = plot_activation_data, annot = annot)

stophere
TESTDIR = 'tests/comparison_prereg'
TESTNAME = 'TRAINING_p_fdr'
THR = 0
myplots.surf.training = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                         plot_function = plot_activation_data)

stophere

