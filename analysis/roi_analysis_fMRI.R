

# load model
rm(list = ls())
cat("\014")

# plot data

# run only for Classic within mask
# plot t1 from where there is an effect
# write how many observations

setwd("~/Software/LeftHand/analysis")
source("./plot_funcs.R")

NCOL = 3

FIGS_DIR = '~/Data/LeftHand/Lund1/figs'
dir.create(FIGS_DIR)
#unlink(file.path(FIGS_DIR, 'figs','*.png'))

THR = 0.99

#################################################################################
radius = 3

surface.rh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/rh.pial.shape.gii'
surface.lh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/lh.pial.shape.gii'
SURFACE_FILES = c(surface.lh, surface.rh)
names(SURFACE_FILES) = c("lh", "rh")

inflated.rh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/rh.inflated.shape.gii'
inflated.lh = '/data/lv0/MotorSkill/labels/subject/fsaverage6/surf/lh.inflated.shape.gii'

DISTANCE = 20

WD = '/data/lv0/MotorSkill/fmriprep/analysis/higherlevel/TrainedCorrect_TrainedIncorrect/'

DATADIR = file.path(WD, 'volume/')

MASK_NAME = 'mask.nii.gz'
TESTDIR = 'tests/linear'
TESTNAME = 'GROUPExp_x_CONDITIONUnt_x_TRAINING_p' #'TRAINING_p_fdr'
myplots.vol.training = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                       plot_function = plot_activation_data)
MASK_NAME = 'mask.nii.gz'
TESTDIR = 'tests/linear'
TESTNAME = 'CONDITIONUnt_p' #'TRAINING_p_fdr'
myplots.vol.training = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                       plot_function = plot_activation_data)

stophere
WD = '/data/lv0/MotorSkill/fmriprep/analysis/higherlevel/Trained_Untrained/'

DATADIR = file.path(WD, 'volume/')

MASK_NAME = 'mask.nii.gz'
TESTDIR = 'tests/linear'
TESTNAME = 'GROUPExp_x_CONDITIONUnt_x_TRAINING_p' #'TRAINING_p_fdr'
myplots.vol.training = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                    plot_function = plot_activation_data)

stophere
radius = 5

DATADIR = file.path(WD, 'surf/')

TESTDIR = 'tests/quadratic'
TESTNAME = 'TRAINING_p_fdr'
THR = 0.95
#myplots.surf.training = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
#                                       plot_function = plot_activation_data)


TESTDIR = 'tests/quadratic'
TESTNAME = 'INTERCEPT_p_fdr'
THR = 0.95
myplots.surf.intercept = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, 
                                         plot_function = plot_activation_data)

