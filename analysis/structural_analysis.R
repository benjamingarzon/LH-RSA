rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)

# study outliers
# covariates

source('~/Software/ImageLMMR/ImageLMMR.R')
setwd("~/Software/LeftHand/analysis")
source('./mytests.R')

SUBJECTS_DIR= '~/Data/LeftHand/Lund1/freesurfer/'
# check 1104
NPROCS = 30

source('./load_covariates.R')

###############################################################
# Define the pipeline
###############################################################

doit = function(WD, MYTEST, OD, 
                MASK = 'rh.thickness.mask.nii.gz', 
                IMAGES_LIST = 'rh.thickness.txt',
                IMAGING_NAME = 'rh.thickness.nii.gz',
                to_gifti = '', 
                NPERMS = 0, 
                shuffle_by = NULL)
  {
  
  setwd(WD)
  IMAGES = read.table(file.path(WD, IMAGES_LIST))
  IMAGING_FILE = file.path(WD, IMAGING_NAME)
  MASK_FILE = file.path(WD, MASK)
  
  OUTPUT_DIR = file.path(WD, OD)
  
  DATA = data.frame(
    SUBJECT = IMAGES$V1, 
    TP = IMAGES$V2
  )
  
  DATA = merge(DATA, covars.table, by = c("SUBJECT", "TP"), all.x = T)
  if (!is.null(IMAGES$V3)) DATA$DEPTH = 1 - IMAGES$V3 # depth starting from cortical surface 

  View(DATA)  
  if ( NPERMS == 0)
  results = vbanalysis(
    IMAGING_FILE,
    OUTPUT_DIR, 
    DATA, 
    MASK_FILE,
    MYTEST,
    remove_outliers = F, to_gifti = to_gifti 
  )
  else 
  results = vbanalysis_perm(
      IMAGING_FILE,
      OUTPUT_DIR, 
      DATA, 
      MASK_FILE,
      MYTEST,
      remove_outliers = F, 
      to_gifti = to_gifti,
      NPERMS = NPERMS,
      shuffle_by = shuffle_by
    )
  
  results$data = DATA
  return(results)
}

###############################################################
# Do it
###############################################################
# get MNI masks for volume tests
# tests T1 x 3 depths, put together as 1 test, 2 hemispheres, quadratic model and linear vs quadratic test
# tests thickness, 2 hemispheres, quadratic model and linear vs quadratic test
# cat 12 volume, quadratic model and linear vs quadratic test
# tests fMRI volume, quadratic model and linear vs quadratic test
# tests fMRI, 2 hemispheres, quadratic model and linear vs quadratic test
# tests FMRI multivariate

#####################
# T1 values
#####################

DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'

rh.gii = file.path(DATADIR, 'rh.fsaverage.white.gii')
lh.gii = file.path(DATADIR, 'lh.fsaverage.white.gii')

rh.mask = '/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/rh.mask.nii.gz'
lh.mask = '/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/lh.mask.nii.gz'
rh.mask = 'rh.mask.nii.gz'
lh.mask = 'lh.mask.nii.gz'

NPERMS = 0

results.reasoning.rh = doit(DATADIR,
                               testreasoning,
                               'tests/thickness/reasoning/',
                               MASK = rh.mask,
                               IMAGES_LIST = 'rh.thickness.txt',
                               IMAGING_NAME = 'rh.thickness.nii.gz',
                               to_gifti = rh.gii,
                               NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))

stophere
results.T1.quadratic.rh = doit(DATADIR,
                               testquadratic_depth,
                               'tests/T1/quadratic.rh',
                               MASK = rh.mask,
                               IMAGES_LIST = 'rh.T1.txt',
                               IMAGING_NAME = 'rh.T1.nii.gz',
                               to_gifti = rh.gii,
                               NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))

results.T1.modelcomparison.rh = doit(DATADIR,
                                    modelcomparison_depth,
                                    'tests/T1/comparison.rh',
                                    MASK = rh.mask,
                                    IMAGES_LIST = 'rh.T1.txt',
                                    IMAGING_NAME = 'rh.T1.nii.gz',
                                    to_gifti = rh.gii,
                                    NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))

results.thickness.quadratic.rh = doit(DATADIR,
                               testquadratic,
                               'tests/thickness/quadratic.rh',
                               MASK = rh.mask,
                               IMAGES_LIST = 'rh.thickness.txt',
                               IMAGING_NAME = 'rh.thickness.nii.gz',
                               to_gifti = rh.gii,
                               NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))

results.thickness.comparison.rh = doit(DATADIR,
                                      modelcomparison,
                                      'tests/thickness/comparison.rh',
                                      MASK = rh.mask,
                                      IMAGES_LIST = 'rh.thickness.txt',
                                      IMAGING_NAME = 'rh.thickness.nii.gz',
                                      to_gifti = rh.gii,
                                      NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))


#################################################
stophere

clusters.stat = read.table('/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results/tests/T1/quadratic-0.75.rh/clusters.stat', header = T)
colnames(clusters.stat) = c("INTERCEPT_p","GROUP_p","TRAINING_p","TRAINING.Q_p",
                            "GROUP_x_TRAINING_p","GROUP_x_TRAINING.Q_p","GROUP_x_TRAINING+_p","GROUP_x_TRAINING-_p",
                            "GROUP_x_TRAINING.Q+_p","GROUP_x_TRAINING.Q-_p")

par(mfrow= c(2, 5))
for(l in colnames(clusters.stat)){
  x = clusters.stat[[l]]
  hist(x, 50, main = l, xlab = "coef") 
  abline(v = x[1], col = 'red')
}
print(clusters.stat)
#####################
# cat12
#####################

DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/cat12'

results.quadratic.cat = doit(DATADIR, 
                                   testquadratic, 
                                   'tests/quadratic', 
                                   MASK = 'mask.nii.gz',  
                                   IMAGES_LIST = 'cat.txt',
                                   IMAGING_NAME = 'data_s8.nii.gz',
                                   NPERMS = 200, shuffle_by = c('SUBJECT', 'TIME')) #, 'INTERCEPT'

stophere



results.quadratic.cat = doit(DATADIR, 
                            testquadratic, 
                            'tests/quadratic.cat', 
                            MASK = 'mask.nii.gz',  
                            IMAGES_LIST = 'cat.txt',
                            IMAGING_NAME = 'data_s8.nii.gz')

results.grouplinearhalf.cat = doit(DATADIR, 
                         testgrouplinearhalf, 
                         'tests/grouplinearhalf', 
                         MASK = 'mask.nii.gz',  
                         IMAGES_LIST = 'cat.txt',
                         IMAGING_NAME = 'data_s8.nii.gz')


DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/vbm/stats'

#####################
# vbm
#####################
# results.grouplinearhalf.vbm = doit(DATADIR, 
#                             testgrouplinearhalf, 
#                             'tests/grouplinearhalf', 
#                             MASK = 'GM_mask.nii.gz',  
#                             IMAGES_LIST = 'vbm.txt',
#                             IMAGING_NAME = 'GM_mod_merg_s4.nii.gz')

results.levels.vbm = doit(DATADIR, 
                                   test2levels, 
                                   'tests/levels', 
                                   MASK = 'GM_mask.nii.gz',  
                                   IMAGES_LIST = 'vbm.txt',
                                   IMAGING_NAME = 'GM_mod_merg_s4.nii.gz')


#####################
# cortical thickness
#####################
DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/'

rh.gii = file.path(DATADIR, 'rh.fsaverage.white.gii')
lh.gii = file.path(DATADIR, 'lh.fsaverage.white.gii')


results.levels.rh = doit(DATADIR,
                                  test2levels,
                                  'tests/levels.rh',
                                  MASK = 'rh.thickness.mask.nii.gz',
                                  IMAGES_LIST = 'rh.thickness.txt',
                                  IMAGING_NAME = 'rh.thickness.nii.gz',
                                  to_gifti = rh.gii)

results.levels.lh = doit(DATADIR,
                         test2levels,
                         'tests/levels.lh',
                         MASK = 'lh.thickness.mask.nii.gz',
                         IMAGES_LIST = 'lh.thickness.txt',
                         IMAGING_NAME = 'lh.thickness.nii.gz',
                         to_gifti = lh.gii)


stophere

results.grouplinearhalf.rh = doit(DATADIR,
                                  testgrouplinearhalf,
                                  'tests/grouplinearhalf.rh',
                                  MASK = 'rh.thickness.mask.nii.gz',
                                  IMAGES_LIST = 'rh.thickness.txt',
                                  IMAGING_NAME = 'rh.thickness.nii.gz',
                                  to_gifti = rh.gii)

results.grouplinearhalf.lh = doit(DATADIR,
                                  testgrouplinearhalf,
                                  'tests/grouplinearhalf.lh',
                                  MASK = 'lh.thickness.mask.nii.gz',
                                  IMAGES_LIST = 'lh.thickness.txt',
                                  IMAGING_NAME = 'lh.thickness.nii.gz',
                                  to_gifti = lh.gii)



stophere
results.quadratic.lh = doit(DATADIR, 
                                  testquadratic, 
                                  'tests/quadratic.lh', 
                                  MASK = 'lh.thickness.mask.nii.gz',  
                                  IMAGES_LIST = 'lh.thickness.txt',
                                  IMAGING_NAME = 'lh.thickness.nii.gz',
                                  to_gifti = lh.gii)

results.quadratic.rh = doit(DATADIR, 
                                  testquadratic, 
                                  'tests/quadratic.rh', 
                                  MASK = 'rh.thickness.mask.nii.gz',  
                                  IMAGES_LIST = 'rh.thickness.txt',
                                  IMAGING_NAME = 'rh.thickness.nii.gz',
                                  to_gifti = rh.gii)

results.anova.lh = doit(DATADIR, 
                            testanova, 
                            'tests/anova.lh', 
                            MASK = 'lh.thickness.mask.nii.gz',  
                            IMAGES_LIST = 'lh.thickness.txt',
                            IMAGING_NAME = 'lh.thickness.nii.gz',
                            to_gifti = lh.gii)

results.anova.rh = doit(DATADIR, 
                            testanova, 
                            'tests/anova.rh', 
                            MASK = 'rh.thickness.mask.nii.gz',  
                            IMAGES_LIST = 'rh.thickness.txt',
                            IMAGING_NAME = 'rh.thickness.nii.gz',
                            to_gifti = rh.gii)

# stophere 
# results.grouplinearhalf.lh = doit(DATADIR, 
#                                   testgrouplinearhalf, 
#                                   'tests/grouplinearhalf.lh', 
#                                   MASK = 'lh.thickness.small.mask.nii.gz',  
#                                   IMAGES_LIST = 'lh.thickness.txt',
#                                   IMAGING_NAME = 'lh.thickness.nii.gz',
#                                   to_gifti = lh.gii)
# 
# results.grouplinearhalf.rh = doit(DATADIR, 
#                                   testgrouplinearhalf, 
#                                   'tests/grouplinearhalf.rh', 
#                                   MASK = 'rh.thickness.small.mask.nii.gz',  
#                                   IMAGES_LIST = 'rh.thickness.txt',
#                                   IMAGING_NAME = 'rh.thickness.nii.gz',
#                                   to_gifti = rh.gii)

vertex = 83215
vertex = 19702
vertex = 151894
#vertex = 90503 # visual
vertex = 139317 # visual
#results = results.quadratic.rh
#results = results.grouplinearhalf.rh
results = results.levels.rh
y = results$imaging.mat[,vertex + 1]
X = results$data
X$y = y
#X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
#X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 

model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
model = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
summary(model)


myplot = ggplot(X, aes(x = TRAINING, group = SUBJECT, col = GROUP, y = y)) + geom_line() + ylim(c(0, 4))
print(myplot)

myplot = ggplot(X, aes(x = TRAINING, group = GROUP, col = GROUP, y = y)) + geom_smooth() + ylim(c(0, 4))
print(myplot)

myplot = ggplot(subset(X, !is.na(TRAINING.L)), aes(x = TRAINING.L, group = GROUP, col = GROUP, y = y)) + geom_point()
print(myplot)

model = lmer(y ~ 1 + GROUP*TRAINING.L + (1 + TRAINING.L|SUBJECT), data = X)
summary(model)

1 - testquadratic(X$y, X)
testanova(X$y, X)

#/home/share/Software/HCP/workbench/bin_rh_linux64/../exe_rh_linux64/wb_command -metric-convert -from-nifti GROUP_x_TRAINING-_p_fdr.nii.gz /home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer//fsaverage/surf/rh.white GROUP_x_TRAINING-_p_fdr.func.gii
stophere


find_outliers = function(X){
  print("Finding outliers")
  mycols = colnames(X)[(grep("\\bV", colnames(X)))]
  l = length(mycols)

  for (k in seq(l)){
    for (subject in unique(X$subject)){
      m = mean(X[X$subject == subject, mycols[k]])
      s = sd(X[X$subject == subject, mycols[k]])
      local_outlier = (abs(X[X$subject == subject, mycols[k]]- m) > 3*s)
      X[local_outlier, mycols[k]] = NA
    }
  }
  
  return(X)
}


########### stophere
results.T1_0.25.quadratic.rh = doit(DATADIR,
                                    testquadratic,
                                    'tests/T1/quadratic-0.25.rh',
                                    MASK = rh.mask,
                                    IMAGES_LIST = 'rh.T1-0.25.txt',
                                    IMAGING_NAME = 'rh.T1-0.25.nii.gz',
                                    to_gifti = rh.gii,
                                    NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))

results.T1_0.75.quadratic.rh = doit(DATADIR,
                                    testquadratic,
                                    'tests/T1/quadratic-0.75.rh',
                                    MASK = rh.mask,
                                    IMAGES_LIST = 'rh.T1-0.75.txt',
                                    IMAGING_NAME = 'rh.T1-0.75.nii.gz',
                                    to_gifti = rh.gii,
                                    NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))

results.T1_0.25.quadratic.lh = doit(DATADIR,
                                    testquadratic,
                                    'tests/T1/quadratic-0.25.lh',
                                    MASK = lh.mask,
                                    IMAGES_LIST = 'lh.T1-0.25.txt',
                                    IMAGING_NAME = 'lh.T1-0.25.nii.gz',
                                    to_gifti = lh.gii,
                                    NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))


results.T1_0.75.quadratic.lh = doit(DATADIR,
                                    testquadratic,
                                    'tests/T1/quadratic-0.75.lh',
                                    MASK = lh.mask,
                                    IMAGES_LIST = 'lh.T1-0.75.txt',
                                    IMAGING_NAME = 'lh.T1-0.75.nii.gz',
                                    to_gifti = lh.gii,
                                    NPERMS = NPERMS, shuffle_by = c('BETWEEN', 'WITHIN'))
