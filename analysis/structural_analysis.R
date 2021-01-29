rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(pracma)
cat("\014")

# check in case there are no clusters
# write out how many observations in each vertex
# do not remove outliers??
# permutation tests
# singular fits and check models
# smoothing
# cat12 w bias
# SPM12
# covariates

# make plots
source('~/Software/ImageLMMR/ImageLMMR.R')
setwd("~/Software/LeftHand/analysis")
source('./mytests.R')

SUBJECTS_DIR= '~/Data/LeftHand/Lund1/freesurfer/'

NPROCS = 35

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
                alpha = 0.05,
                shuffle_by = NULL, upsample = NULL, remove_outliers = T)
  {
  
  setwd(WD)
  IMAGES = read.table(file.path(WD, IMAGES_LIST))
  IMAGING_FILE = file.path(WD, IMAGING_NAME)
  MASK_FILE = file.path(WD, MASK)
  
  OUTPUT_DIR = file.path(WD, OD)

  if (is.null(IMAGES$V3)) {
    DATA = data.frame(
      SUBJECT = IMAGES$V1, 
      TP = IMAGES$V2
    )
    DATA = merge(DATA, covars.table, by = c("SUBJECT", "TP"), all.x = T) %>% arrange(SUBJECT, TP)
    
  } 
  else {
    DATA = data.frame(
      SUBJECT = IMAGES$V1, 
      TP = IMAGES$V2,
      DEPTH = 1 - IMAGES$V3 # depth starting from cortical surface/ inverse depth
    )
    DATA = merge(DATA, covars.table, by = c("SUBJECT", "TP"), all.x = T) %>% arrange(-DEPTH, SUBJECT, TP)
    
  }

  View(DATA)
#  browser()
  excluded = which((DATA$SUBJECT %in% names(which(table(DATA$SUBJECT) < 4))) | !complete.cases(DATA))
  if ( NPERMS == 0)
  results = vbanalysis(
    IMAGING_FILE,
    OUTPUT_DIR, 
    DATA, 
    MASK_FILE,
    MYTEST,
    remove_outliers = remove_outliers, 
    excluded = excluded, 
    to_gifti = to_gifti, alpha = alpha, upsample = upsample,
  )
  else 
  results = vbanalysis_perm(
      IMAGING_FILE,
      OUTPUT_DIR, 
      DATA, 
      MASK_FILE,
      MYTEST,
      remove_outliers = remove_outliers, 
      excluded = excluded, 
      to_gifti = to_gifti, alpha = alpha, upsample = upsample,
      NPERMS = NPERMS,
      shuffle_by = shuffle_by
    )
  
  results$data = DATA
  save(results, file = file.path(OUTPUT_DIR, 'results.rda'))
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



FSDIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'
CATDIR = '/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces/' #_32k/' 
rh.gii = file.path(FSDIR, 'rh.fsaverage.white.gii')
lh.gii = file.path(FSDIR, 'lh.fsaverage.white.gii')
rh.cat.gii = file.path(CATDIR, 'rh.central.freesurfer.gii')
lh.cat.gii = file.path(CATDIR, 'lh.central.freesurfer.gii')

#rh.mask = '/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/rh.mask.nii.gz'
#lh.mask = '/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage/lh.mask.nii.gz'
rh.mask = 'rh.mask.nii.gz'
lh.mask = 'lh.mask.nii.gz'
rh.cortex.mask = 'rh.cortex.mask.nii.gz'
lh.cortex.mask = 'lh.cortex.mask.nii.gz'

#vbm.mask = 'mask_whole_brain_3mm.nii.gz'
#vbm.data = 'data_s8_3mm.nii.gz'
vbm.mask = 'mask.nii.gz' 
#vbm.mask = 'mask_whole_brain.nii.gz'
vbm.data = 'data_s8.nii.gz'

upsample = NULL #'/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/mask_whole_brain.nii.gz'

NPERMS = 50
shuffle_by = c('BETWEEN', 'WITHIN')

alphavoxel = 0.05
alphavertex = 0.025

#####################
# cat12 surf
#####################

DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/cat12crossbias8_10'
if (T){
results.thickness.cat.quadratic.rh = doit(DATADIR,
                                           testquadratic,
                                           'tests/quadratic.rh',
                                           MASK = rh.mask,
                                           IMAGES_LIST = 'rh.thickness.txt',
                                           IMAGING_NAME = 'rh.thickness.10.nii.gz',
                                           to_gifti = rh.cat.gii, alpha = alphavertex,
                                           NPERMS = 0, shuffle_by = shuffle_by)

results.thickness.cat.quadratic.lh = doit(DATADIR,
                                           testquadratic,
                                           'tests/quadratic.lh',
                                           MASK = lh.mask,
                                           IMAGES_LIST = 'lh.thickness.txt',
                                           IMAGING_NAME = 'lh.thickness.10.nii.gz',
                                           to_gifti = lh.cat.gii, alpha = alphavertex,
                                           NPERMS = 0, shuffle_by = shuffle_by)


results.thickness.cat.comparison.rh = doit(DATADIR,
                                           modelcomparison,
                                           'tests/comparison.rh',
                                           MASK = rh.mask,
                                           IMAGES_LIST = 'rh.thickness.txt',
                                           IMAGING_NAME = 'rh.thickness.10.nii.gz',
                                           to_gifti = rh.cat.gii, alpha = alphavertex,
                                           NPERMS = 0, shuffle_by = shuffle_by)


results.thickness.cat.comparison.lh = doit(DATADIR,
                                           modelcomparison,
                                           'tests/comparison.lh',
                                           MASK = lh.mask,
                                           IMAGES_LIST = 'lh.thickness.txt',
                                           IMAGING_NAME = 'lh.thickness.10.nii.gz',
                                           to_gifti = lh.cat.gii, alpha = alphavertex,
                                           NPERMS = 0, shuffle_by = shuffle_by)

results.thickness.cat.reliability.rh = doit(DATADIR,
                                           reliability,
                                           'tests/reliability.rh',
                                           MASK = rh.cortex.mask,
                                           IMAGES_LIST = 'rh.thickness.txt',
                                           IMAGING_NAME = 'rh.thickness.10.nii.gz',
                                           to_gifti = rh.cat.gii, alpha = alphavertex,
                                           NPERMS = 0, shuffle_by = shuffle_by)

}
#####################
# cat12
#####################

#DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/cat12crossbias12_15'
DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/cat12crossbias8_10'

results.reliability.cat = doit(DATADIR, 
                             reliability, 
                             'tests/reliability', 
                             MASK = 'mask_whole_brain.nii.gz',
                             IMAGES_LIST = 'image_list.txt',
                             IMAGING_NAME = vbm.data, upsample = upsample, 
                             NPERMS = 0, shuffle_by = shuffle_by) 
stophere
results.quadratic.cat = doit(DATADIR, 
                             testquadratic, 
                             'tests/quadratic', 
                             MASK = vbm.mask,  
                             IMAGES_LIST = 'image_list.txt',
                             IMAGING_NAME = vbm.data, upsample = upsample, 
                             NPERMS = NPERMS, shuffle_by = shuffle_by) 

results.comparison.cat = doit(DATADIR, 
                             modelcomparison, 
                             'tests/comparison', 
                             MASK = vbm.mask,  
                             IMAGES_LIST = 'image_list.txt',
                             IMAGING_NAME = vbm.data, upsample = upsample, 
                             NPERMS = 0, shuffle_by = shuffle_by) 



if (F){
  
results.linear.cat = doit(DATADIR,
                              testlinear,
                              'tests/linear',
                              MASK = vbm.mask,
                              IMAGES_LIST = 'image_list.txt',
                              IMAGING_NAME = vbm.data, upsample = upsample,
                              NPERMS = 0, shuffle_by = shuffle_by)


results.reasoning.cat = doit(DATADIR,
                             testreasoning,
                             'tests/reasoning',
                             MASK = vbm.mask,  
                             IMAGES_LIST = 'image_list.txt',
                             IMAGING_NAME = vbm.data, upsample = upsample, 
                             NPERMS = 0, shuffle_by = shuffle_by) 

}

#####################
# thickness
#####################

DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'

if (F){

  results.thickness.reliability.rh = doit(DATADIR,
                                          reliability,
                                          'tests/thickness/reliability.rh',
                                          MASK = rh.cortex.mask,
                                          IMAGES_LIST = 'rh.thickness.txt',
                                          IMAGING_NAME = 'rh.thickness.nii.gz',
                                          to_gifti = rh.gii, alpha = alphavertex,
                                          NPERMS = NPERMS, shuffle_by = shuffle_by, remove_outliers = F)
  
  results.thickness.reliability.lh = doit(DATADIR,
                                          reliability,
                                          'tests/thickness/reliability.lh',
                                          MASK = lh.cortex.mask,
                                          IMAGES_LIST = 'lh.thickness.txt',
                                          IMAGING_NAME = 'lh.thickness.nii.gz',
                                          to_gifti = lh.gii, alpha = alphavertex,
                                          NPERMS = NPERMS, shuffle_by = shuffle_by, remove_outliers = F)
  
results.thickness.linear.rh = doit(DATADIR,
                                       testlinear,
                                       'tests/thickness/linear.rh',
                                       MASK = rh.mask,
                                       IMAGES_LIST = 'rh.thickness.txt',
                                       IMAGING_NAME = 'rh.thickness.nii.gz',
                                       to_gifti = rh.gii, alpha = alphavertex,
                                       NPERMS = NPERMS, shuffle_by = shuffle_by)

results.thickness.linear.lh = doit(DATADIR,
                                       testlinear,
                                       'tests/thickness/linear.lh',
                                       MASK = lh.mask,
                                       IMAGES_LIST = 'lh.thickness.txt',
                                       IMAGING_NAME = 'lh.thickness.nii.gz',
                                       to_gifti = lh.gii, alpha = alphavertex,
                                       NPERMS = NPERMS, shuffle_by = shuffle_by)
  

}

#right
results.thickness.quadratic.rh = doit(DATADIR,
                                      testquadratic, #_Classic,
                                      'tests/thickness/quadratic.rh',
                                      MASK = rh.mask,
                                      IMAGES_LIST = 'rh.thickness.txt',
                                      IMAGING_NAME = 'rh.thickness.nii.gz',
                                      to_gifti = rh.gii, alpha = alphavertex, 
                                      NPERMS = NPERMS, shuffle_by = shuffle_by)

# left
results.thickness.quadratic.lh = doit(DATADIR,
                                      testquadratic, #_Classic,
                                      'tests/thickness/quadratic.lh',
                                      MASK = lh.mask,
                                      IMAGES_LIST = 'lh.thickness.txt',
                                      IMAGING_NAME = 'lh.thickness.nii.gz',
                                      to_gifti = lh.gii, alpha = alphavertex, 
                                      NPERMS = NPERMS, shuffle_by = shuffle_by)


results.thickness.comparison.rh = doit(DATADIR,
                                       modelcomparison,
                                       'tests/thickness/comparison.rh',
                                       MASK = rh.mask,
                                       IMAGES_LIST = 'rh.thickness.txt',
                                       IMAGING_NAME = 'rh.thickness.nii.gz',
                                       to_gifti = rh.gii, alpha = alphavertex,
                                       NPERMS = 0, shuffle_by = shuffle_by)


results.thickness.comparison.lh = doit(DATADIR,
                                       modelcomparison,
                                       'tests/thickness/comparison.lh',
                                       MASK = rh.mask,
                                       IMAGES_LIST = 'lh.thickness.txt',
                                       IMAGING_NAME = 'lh.thickness.nii.gz',
                                       to_gifti = lh.gii, alpha = alphavertex,
                                       NPERMS = 0, shuffle_by = shuffle_by)




#####################
# T1 values
#####################

# right
results.T1.quadratic.rh = doit(DATADIR,
                               testquadratic_depth,
                               'tests/T1/quadratic.rh',
                               MASK = rh.mask,
                               IMAGES_LIST = 'rh.T1.txt',
                               IMAGING_NAME = 'rh.T1.nii.gz',
                               to_gifti = rh.gii,  alpha = alphavertex,
                               NPERMS = NPERMS, shuffle_by = shuffle_by)

#left
results.T1.quadratic.lh = doit(DATADIR,
                               testquadratic_depth,
                               'tests/T1/quadratic.lh',
                               MASK = lh.mask,
                               IMAGES_LIST = 'lh.T1.txt',
                               IMAGING_NAME = 'lh.T1.nii.gz',
                               to_gifti = lh.gii, alpha = alphavertex,
                               NPERMS = NPERMS, shuffle_by = shuffle_by)

if (F){
  
results.T1.linearhalf.rh = doit(DATADIR,
                               testgrouplinearhalf_depth,
                               'tests/T1/linearhalf.rh',
                               MASK = rh.mask,
                               IMAGES_LIST = 'rh.T1.txt',
                               IMAGING_NAME = 'rh.T1.nii.gz',
                               to_gifti = rh.gii,  alpha = alphavertex,
                               NPERMS = NPERMS, shuffle_by = shuffle_by)

results.T1.linearhalf.lh = doit(DATADIR,
                            testgrouplinearhalf_depth,
                            'tests/T1/linearhalf.lh',
                            MASK = lh.mask,
                            IMAGES_LIST = 'lh.T1.txt',
                            IMAGING_NAME = 'lh.T1.nii.gz',
                            to_gifti = lh.gii,  alpha = alphavertex,
                            NPERMS = NPERMS, shuffle_by = shuffle_by)

  

results.T1.modelcomparison.rh = doit(DATADIR,
                                    modelcomparison_depth,
                                    'tests/T1/comparison.rh',
                                    MASK = rh.mask,
                                    IMAGES_LIST = 'rh.T1.txt',
                                    IMAGING_NAME = 'rh.T1.nii.gz',
                                    to_gifti = rh.gii,  alpha = alphavertex, 
                                    NPERMS = NPERMS, shuffle_by = shuffle_by)

results.T1.modelcomparison.lh = doit(DATADIR,
                                     modelcomparison_depth,
                                     'tests/T1/comparison.lh',
                                     MASK = lh.mask,
                                     IMAGES_LIST = 'lh.T1.txt',
                                     IMAGING_NAME = 'lh.T1.nii.gz',
                                     to_gifti = lh.gii, alpha = alphavertex, 
                                     NPERMS = NPERMS, shuffle_by = shuffle_by)

}

stophere

results.thickness.variance.rh = doit(DATADIR,
                                     testvariance,
                                     'tests/thickness/variance.rh',
                                     MASK = rh.mask,
                                     IMAGES_LIST = 'rh.thickness.txt',
                                     IMAGING_NAME = 'rh.thickness.nii.gz',
                                     to_gifti = rh.gii, alpha = alphavertex, 
                                     NPERMS = NPERMS, shuffle_by = shuffle_by)

results.thickness.variance.lh = doit(DATADIR,
                                     testvariance,
                                     'tests/thickness/variance.lh',
                                     MASK = lh.mask,
                                     IMAGES_LIST = 'lh.thickness.txt',
                                     IMAGING_NAME = 'lh.thickness.nii.gz',
                                     to_gifti = lh.gii, alpha = alphavertex, 
                                     NPERMS = NPERMS, shuffle_by = shuffle_by)

#####################
# GMV
#####################


DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/cat12'

results.variance.cat = doit(DATADIR, 
                            testvariance, 
                            'tests/variance', 
                            MASK = vbm.mask,  
                            IMAGES_LIST = 'image_list.txt',
                            IMAGING_NAME = vbm.data,
                            NPERMS = 0, shuffle_by = shuffle_by) 

#####################
# Other tests
#####################

if (F)
  #####################
# Reasoning
#####################
DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'

# right
results.reasoning.rh = doit(DATADIR,
                            testreasoning,
                            'tests/thickness/reasoning.rh',
                            MASK = rh.mask,
                            IMAGES_LIST = 'rh.thickness.txt',
                            IMAGING_NAME = 'rh.thickness.nii.gz',
                            to_gifti = rh.gii,  alpha = alphavertex, 
                            NPERMS = NPERMS, shuffle_by = shuffle_by)

# left
results.reasoning.lh = doit(DATADIR,
                            testreasoning,
                            'tests/thickness/reasoning.lh',
                            MASK = lh.mask,
                            IMAGES_LIST = 'lh.thickness.txt',
                            IMAGING_NAME = 'lh.thickness.nii.gz',
                            to_gifti = lh.gii,  alpha = alphavertex, 
                            NPERMS = NPERMS, shuffle_by = shuffle_by)


#####################
# Coil
#####################








stophere



results.quadratic.cat = doit(DATADIR, 
                            testquadratic, 
                            'tests/quadratic.cat', 
                            MASK = vbm.mask,  
                            IMAGES_LIST = 'cat.txt',
                            IMAGING_NAME = vbm.data)

results.grouplinearhalf.cat = doit(DATADIR, 
                         testgrouplinearhalf, 
                         'tests/grouplinearhalf', 
                         MASK = vbm.mask,  
                         IMAGES_LIST = 'cat.txt',
                         IMAGING_NAME = vbm.data)


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

vertex = 832
vertex = which.max(results$stats[, 'GROUP_x_TRAINING_p'])
vertex = which.min(abs(results$stats[, 'GROUP_x_TRAINING_tstat'] + 3.072283506))
results = results.thickness.variance.rh # results.linear.cat

y = results$imaging.mat[,vertex + 1]
X = results$data
X$y = y
X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
#X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 

model = lmer(y ~ 1 + SYSTEM + GROUP*TRAINING + (1 + TRAINING|SUBJECT), data = X)

#model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
#model = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
summary(model)


myplot = ggplot(X, aes(x = TRAINING, group = SUBJECT, col = GROUP, y = y)) + geom_line() + geom_point() + ylim(c(0, 1))
print(myplot)

myplot = ggplot(X, aes(x = TRAINING, group = GROUP, col = GROUP, y = y)) + geom_smooth() + ylim(c(0.2, .7))
print(myplot)

myplot = ggplot(subset(X, !is.na(TRAINING.L)), aes(x = TRAINING.L, group = GROUP, col = GROUP, y = y)) + geom_point()
print(myplot)

model = lmer(y ~ 1 + GROUP*TRAINING.L + (1 + TRAINING.L|SUBJECT), data = X)
summary(model)

1 - testquadratic(X$y, X)
testanova(X$y, X)

#/home/share/Software/HCP/workbench/bin_rh_linux64/../exe_rh_linux64/wb_command -metric-convert -from-nifti GROUP_x_TRAINING-_p_fdr.nii.gz /home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer//fsaverage/surf/rh.white GROUP_x_TRAINING-_p_fdr.func.gii
stophere


# find_outliers = function(X){
#   print("Finding outliers")
#   mycols = colnames(X)[(grep("\\bV", colnames(X)))]
#   l = length(mycols)
# 
#   for (k in seq(l)){
#     for (subject in unique(X$subject)){
#       m = mean(X[X$subject == subject, mycols[k]])
#       s = sd(X[X$subject == subject, mycols[k]])
#       local_outlier = (abs(X[X$subject == subject, mycols[k]]- m) > 3*s)
#       X[local_outlier, mycols[k]] = NA
#     }
#   }
#   
#   return(X)
# }


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
