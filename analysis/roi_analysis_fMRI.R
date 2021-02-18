

# load model
rm(list = ls())
cat("\014")

# plot data

# run only for Classic within mask
# plot t1 from where there is an effect
# write how many observations

setwd("~/Software/LeftHand/analysis")
source("./plot_funcs.R")

NCOL = 4

FIGS_DIR = '~/Data/LeftHand/Lund1/figs'
dir.create(FIGS_DIR)
#unlink(file.path(FIGS_DIR, 'figs','*.png'))

THR = 0.95



#################################################################################
radius = 10


surface.rh = '~/Data/LeftHand/Lund1/labels/subject/fsaverage6/surf/rh.pial.gii'
surface.lh = '~/Data/LeftHand/Lund1/labels/subject/fsaverage6/surf/lh.pial.gii'
SURFACE_FILES = c(surface.lh, surface.rh)
names(SURFACE_FILES) = c("lh", "rh")

inflated.rh = '~/Data/LeftHand/Lund1/labels/subject/fsaverage6/surf/rh.inflated'
inflated.lh = '~/Data/LeftHand/Lund1/labels/subject/fsaverage6/surf/lh.inflated'


DISTANCE = 20
DATADIR = '/data/lv0/MotorSkill/fmriprep/analysis/higherlevel/Trained_Untrained/surfR/'

MASK_NAME = 'mask.nii.gz'
# thickness linear
TESTDIR = 'tests/linear'
TESTNAME = 'GROUP_x_TRAINING.Q-_p'
myplots.thickness.linear.groupxtraining = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR)

#stophere
#TESTNAME = 'GROUP_x_TRAINING.Q_p'
#myplots.thickness.quadratic.groupxtraining = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR)

TESTDIR = 'tests/thickness/quadratic'
TESTNAME = 'Omni_p'
#THR = 0.95
#myplots.thickness.quadratic.omni = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR)

TESTDIR = 'tests/T1/quadratic'
TESTNAME = 'thickness_GROUP_x_TRAINING-_p'
myplots.T1.quadratic.groupxtraining = create_surf_rois(
  DATADIR,
  TESTDIR,
  TESTNAME,
  DISTANCE,
  radius,
  MASK_NAME,
  THR,
  myplots.thickness.linear.groupxtraining$ROI_FILE,
  wDEPTH = T
)

