# load model
rm(list = ls())
cat("\014")

# plot data
#ROIS: 
#VBM: 1 hand motor, 2: PCC, 3: visual
#thickness: rh: 1 hand motor, 2: PCC, 3: visual, lh : 1 visual, 2 parietal
# only hand motor coincides for thickness and VBM

setwd("~/Software/LeftHand/analysis")
source("./plot_funcs.R")

NCOL = 6

FIGS_DIR = '~/Data/LeftHand/Lund1/figs/structure'
dir.create(FIGS_DIR)

GMV_LABEL = 'GMV (a. u.)' #% increase in GMV'
THICKNESS_LABEL = 'Cortical thickness (mm)' # % increase in cortical thickness'
#################################################################################
# VBM
radius = 4
DISTANCE = 10
DATADIR = '~/Data/LeftHand/Lund1/cat12/'
MASK_NAME = 'mask.nii.gz'
THR=0.95
REDRAW = T

TESTDIR = 'tests/linear'
TESTNAME = 'GROUP_x_TRAINING_p_fdr'
myplots.VBM = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, GMV_LABEL, redraw = REDRAW)


TESTDIR = 'tests/linear'
TESTNAME = 'GROUP+_p'
myplots.VBM = create_vol_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, GMV_LABEL, redraw = REDRAW)


#################################################################################
# thickness

surface.rh='/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/rh.pial.gii'
surface.lh='/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/lh.pial.gii'

SURFACE_FILES = c(surface.lh, surface.rh)
names(SURFACE_FILES) = c("lh", "rh")

inflated.rh='/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/rh.inflated'
inflated.lh='/home/benjamin.garzon/Data/LeftHand/Lund1/labels/subject/fsaverage/surf/lh.inflated'

DISTANCE=20
radius=5
DATADIR = '/home/share/MotorSkill/freesurfer/results/'

MASK_NAME = 'mask.nii.gz'
THR = 0.975

# thickness linear
TESTDIR = 'tests/thickness/linear'
TESTNAME = 'GROUP_x_TRAINING_p'
#myplots.thickness.GROUP_x_TRAINING_p = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, THICKNESS_LABEL, redraw = REDRAW)

TESTNAME = 'GROUP_p'
#myplots.thickness.GROUP_p = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, THICKNESS_LABEL, redraw = REDRAW)

# T1 asymptotic

DATADIR = '/home/share/MotorSkill/freesurfer/results/'
TESTDIR = 'tests/T1/quadratic'
TESTNAME = 'GROUP_x_TRAINING_p'
myplots.T1.GROUP_x_TRAINING_p = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, THICKNESS_LABEL, redraw = REDRAW)

TESTNAME = 'GROUP_x_TRAINING.Q_p'
myplots.T1.GROUP_x_TRAINING.Q_p = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, THICKNESS_LABEL, redraw = REDRAW)

TESTNAME = 'GROUP_p'
myplots.T1.GROUP_p = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, THICKNESS_LABEL, redraw = REDRAW)

stophere

# thickness comparison
TESTDIR = 'tests/thickness/comparison'
TESTNAME = 'Best'
myplots.thickness.comparison = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, THICKNESS_LABEL, redraw = REDRAW)


stophere
#################################################################################
# cat thickness

surface.rh='/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces_32k/rh.central.freesurfer.gii'
surface.lh='/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces_32k/lh.central.freesurfer.gii'

DISTANCE=20
radius=5
DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/cat12cross/'
FIGS_DIR = file.path(DATADIR, 'figs')
dir.create(FIGS_DIR)
TESTNAME = 'TRAINING+_p.func.gii'

MASK_NAME = 'cortex.mask.nii.gz'
TESTDIR = 'tests/linear'

SURFACE_FILES = c(surface.lh, surface.rh)
names(SURFACE_FILES) = c("lh", "rh")
#myplots = create_surf_rois(DATADIR, TESTDIR, DISTANCE, radius, MASK_NAME)
#print(ggarrange(plotlist = myplots, nrow = ceiling(length(myplots)/NCOL), ncol = NCOL))



stophere

TESTNAME = 'thickness_Omni_p'
myplots.T1.quadratic.omni = create_surf_rois(
  DATADIR,
  TESTDIR,
  TESTNAME,
  DISTANCE,
  radius,
  MASK_NAME,
  THR,
  myplots.thickness.quadratic.omni$ROI_FILE,
  wDEPTH = T
)

TESTDIR = 'tests/T1/quadratic'
TESTNAME = 'GROUP_x_TRAINING_p'
myplots.T1.quadratic.groupxtraining = create_surf_rois(DATADIR,
                                                       TESTDIR,
                                                       TESTNAME,
                                                       DISTANCE,
                                                       radius,
                                                       MASK_NAME,
                                                       THR,
                                                       wDEPTH = T)


#TESTDIR = 'tests/T1/quadratic'
#TESTNAME = 'GROUP_x_TRAINING.Q_p.func.gii'
#myplots.T1.quadratic.groupxtraining = create_surf_rois(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR)

stophere

##############################

mask <- fast_readnii(MASK_FILE)
whichroi = 1

roi <- fast_readnii(ROI_FILE)[, , , whichroi]
roi_indices = which(roi[mask > 0] > 0)

#roi_indices = seq(100)
imaging.mat = results$imaging.mat[, roi_indices]

X = results$data[-results$excluded,]
X$y = rowMeans(imaging.mat)

model.l = lmer(y ~ 1 + SYSTEM + GROUP * TRAINING * DEPTH + (1 |
                                                              SUBJECT), data = X)
model.q = lmer(y ~ 1 + SYSTEM + GROUP * (TRAINING + TRAINING.Q) + (1 |
                                                                     SUBJECT), data = X)

summary(model.l)
summary(model.q)
View(X)

myplot = ggplot(X, aes(
  x = TP,
  group = SUBJECT,
  col = GROUP,
  y = y
)) + geom_line() + geom_point() + facet_grid(. ~ DEPTH) +
  scale_colour_manual(values = myPalette) + theme_lh#+ ylim(c(0, 1))
print(myplot)

myplot = ggplot(subset(X,!is.na(y)), aes(x = GROUP, fill = GROUP)) + geom_bar() + facet_grid(DEPTH ~ TP) +
  scale_fill_manual(values = myPalette) + theme_lh#+ ylim(c(0, 1))
print(myplot)

#myplot = ggplot(X, aes(x = TP, group = GROUP, col = GROUP, y = y)) + geom_smooth() + ylim(c(0.2, .4))
#print(myplot)


X.sem = X %>% group_by(TP, GROUP) %>% dplyr::summarise(y.mean = mean(y, na.rm = T), y.sem = sem(y))
myplot = ggplot(
  X.sem,
  aes(
    x = TP,
    group = GROUP,
    col = GROUP,
    y = y.mean,
    ymax = y.mean + y.sem,
    ymin = y.mean - y.sem
  )
) +
  geom_line() + geom_point() + geom_errorbar() +
  scale_colour_manual(values = myPalette) + theme_lh #+ ylim(c(0.2, .3))
print(myplot)

X.sem = X %>% group_by(TP, GROUP, DEPTH) %>% dplyr::summarise(y.mean = mean(y, na.rm = T), y.sem = sem(y))

myplot = ggplot(
  X.sem,
  aes(
    x = TP,
    group = GROUP,
    col = GROUP,
    y = y.mean,
    ymax = y.mean + y.sem,
    ymin = y.mean - y.sem
  )
) +
  geom_line() + geom_point() + geom_errorbar() + facet_grid(. ~ DEPTH)
print(myplot)
