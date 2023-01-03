rm(list = ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(pracma)
cat("\014")

# do not remove outliers??
# singular fits and check models
# cat12 w bias
# SPM12
# covariates

# make plots
source('~/Software/ImageLMMR/ImageLMMR.R')
setwd("~/Software/LeftHand/analysis")
source('./mytests.R')
source('./structural_analysis_funcs.R')

SUBJECTS_DIR = '~/Data/LeftHand/Lund1/freesurfer/'
NPROCS = 20

source('./load_covariates.R')

###############################################################
# Define the pipeline
###############################################################
doit = structural_analysis_univariate

FSDIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'
CATDIR = '/home/benjamin.garzon/Software/spm/spm12/toolbox/cat12/templates_surfaces/' #_32k/'
rh.gii = file.path(FSDIR, 'rh.fsaverage.white.gii')
lh.gii = file.path(FSDIR, 'lh.fsaverage.white.gii')
rh.cat.gii = file.path(CATDIR, 'rh.central.freesurfer.gii')
lh.cat.gii = file.path(CATDIR, 'lh.central.freesurfer.gii')

rh.mask = 'rh.mask.nii.gz'
lh.mask = 'lh.mask.nii.gz'
rh.cortex.mask = 'rh.cortex.mask.nii.gz'
lh.cortex.mask = 'lh.cortex.mask.nii.gz'

vbm.mask = 'mask.nii.gz' #vbm.mask = 'mask_whole_brain_3mm.nii.gz'
#vbm.mask = 'mask_whole_brain.nii.gz'
vbm.data = 'data_s8.nii.gz' #vbm.data = 'data_s8_3mm.nii.gz'

upsample = NULL #'/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/mask_whole_brain.nii.gz'

# only for perm tests
NPERMS = 0
shuffle_by = c('BETWEEN', 'WITHIN')

alphavoxel = 0.05
alphavertex = 0.025

remove_outliers = F
###############################################################
# Do it
###############################################################

#####################
# cat12
#####################

DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/cat12'

if (F){
#results.reliability.cat =
doit(
  DATADIR,
  reliability,
  'tests/reliability',
  MASK = 'mask_whole_brain.nii.gz',
  IMAGES_LIST = 'image_list.txt',
  IMAGING_NAME = vbm.data,
  upsample = upsample,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

doit(
    DATADIR,
    reliability,
    'tests/reliability-mask',
    MASK = vbm.mask,
    IMAGES_LIST = 'image_list.txt',
    IMAGING_NAME = vbm.data,
    upsample = upsample,
    NPERMS = NPERMS,
    shuffle_by = shuffle_by
  )

#results.quadratic.cat =
doit(
  DATADIR,
  testquadratic,
  'tests/quadratic',
  MASK = vbm.mask,
  IMAGES_LIST = 'image_list.txt',
  IMAGING_NAME = vbm.data,
  upsample = upsample,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.linear.cat =
doit(
  DATADIR,
  testlinear,
  'tests/linear',
  MASK = vbm.mask,
  IMAGES_LIST = 'image_list.txt',
  IMAGING_NAME = vbm.data,
  upsample = upsample,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.asymptotic.cat =
doit(
  DATADIR,
  testasymptotic,
  'tests/asymptotic',
  MASK = vbm.mask,
  IMAGES_LIST = 'image_list.txt',
  IMAGING_NAME = vbm.data,
  upsample = upsample,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.comparison.cat =
doit(
  DATADIR,
  modelcomparison,
  'tests/comparison',
  MASK = vbm.mask,
  IMAGES_LIST = 'image_list.txt',
  IMAGING_NAME = vbm.data,
  upsample = upsample,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

}

#####################
# thickness
#####################

DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'

if (F){
#results.thickness.reliability.rh =
doit(
  DATADIR,
  reliability,
  'tests/thickness/reliability.rh',
  MASK = rh.cortex.mask,
  IMAGES_LIST = 'rh.thickness.txt',
  IMAGING_NAME = 'rh.thickness.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.reliability.lh =
doit(
  DATADIR,
  reliability,
  'tests/thickness/reliability.lh',
  MASK = lh.cortex.mask,
  IMAGES_LIST = 'lh.thickness.txt',
  IMAGING_NAME = 'lh.thickness.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.reliability.rh =
doit(
  DATADIR,
  reliability,
  'tests/thickness/reliability-mask.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.thickness.txt',
  IMAGING_NAME = 'rh.thickness.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.reliability.lh =
doit(
  DATADIR,
  reliability,
  'tests/thickness/reliability-mask.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.thickness.txt',
  IMAGING_NAME = 'lh.thickness.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.quadratic.rh =
doit(
  DATADIR,
  testquadratic,
  'tests/thickness/quadratic.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.thickness.txt',
  IMAGING_NAME = 'rh.thickness.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.quadratic.lh =
doit(
  DATADIR,
  testquadratic,
  'tests/thickness/quadratic.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.thickness.txt',
  IMAGING_NAME = 'lh.thickness.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.linear.rh =
doit(
  DATADIR,
  testlinear,
  'tests/thickness/linear.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.thickness.txt',
  IMAGING_NAME = 'rh.thickness.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.linear.lh =
doit(
  DATADIR,
  testlinear,
  'tests/thickness/linear.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.thickness.txt',
  IMAGING_NAME = 'lh.thickness.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)


#results.thickness.asymptotic.rh =
doit(
  DATADIR,
  testasymptotic,
  'tests/thickness/asymptotic.rh',
  MASK = lh.mask,
  IMAGES_LIST = 'rh.thickness.txt',
  IMAGING_NAME = 'rh.thickness.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.asymptotic.lh =
doit(
  DATADIR,
  testasymptotic,
  'tests/thickness/asymptotic.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.thickness.txt',
  IMAGING_NAME = 'lh.thickness.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)


#results.thickness.comparison.rh =
doit(
  DATADIR,
  modelcomparison,
  'tests/thickness/comparison.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.thickness.txt',
  IMAGING_NAME = 'rh.thickness.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.thickness.comparison.lh =
doit(
  DATADIR,
  modelcomparison,
  'tests/thickness/comparison.lh',
  MASK = rh.mask,
  IMAGES_LIST = 'lh.thickness.txt',
  IMAGING_NAME = 'lh.thickness.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

}
#####################
# T1 values
#####################


#asymptotic
#right
doit(
  DATADIR,
  testasymptotic_depth,
  'tests/T1/asymptotic.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.T1.txt',
  IMAGING_NAME = 'rh.T1.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#left
#results.T1.asymptotic.lh = 
doit(
  DATADIR,
  testasymptotic_depth,
  'tests/T1/asymptotic.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.T1.txt',
  IMAGING_NAME = 'lh.T1.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.T1.quadratic.rh = 
doit(
  DATADIR,
  testcubic_depth,
  'tests/T1/cubic.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.T1.txt',
  IMAGING_NAME = 'rh.T1.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#left
#results.T1.quadratic.lh = 
doit(
  DATADIR,
  testcubic_depth,
  'tests/T1/cubic.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.T1.txt',
  IMAGING_NAME = 'lh.T1.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

stophere


if (F) {

  #model comparison
  #right
  #results.T1.modelcomparison.rh = 
  doit(
    DATADIR,
    modelcomparison_depth,
    'tests/T1/comparison.rh',
    MASK = rh.mask,
    IMAGES_LIST = 'rh.T1.txt',
    IMAGING_NAME = 'rh.T1.nii.gz',
    to_gifti = rh.gii,
    alpha = alphavertex,
    NPERMS = NPERMS,
    shuffle_by = shuffle_by
  )
  
  #left
  #results.T1.modelcomparison.lh = 
  doit(
    DATADIR,
    modelcomparison_depth,
    'tests/T1/comparison.lh',
    MASK = lh.mask,
    IMAGES_LIST = 'lh.T1.txt',
    IMAGING_NAME = 'lh.T1.nii.gz',
    to_gifti = lh.gii,
    alpha = alphavertex,
    NPERMS = NPERMS,
    shuffle_by = shuffle_by
  )
  
  

# linear
#right
#results.T1.linear.rh = 
doit(
  DATADIR,
  testlinear_depth,
  'tests/T1/linear.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.T1.txt',
  IMAGING_NAME = 'rh.T1.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#left
#results.T1.linear.lh = 
doit(
  DATADIR,
  testlinear_depth,
  'tests/T1/linear.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.T1.txt',
  IMAGING_NAME = 'lh.T1.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)


#reliability
#right

#results.T1.reliability.rh =
doit(
  DATADIR,
  reliability_depth,
  'tests/T1/reliability.rh',
  MASK = rh.cortex.mask,
  IMAGES_LIST = 'rh.T1.txt',
  IMAGING_NAME = 'rh.T1.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#results.T1.reliability.lh =
doit(
  DATADIR,
  reliability_depth,
  'tests/T1/reliability-mask.rh',
  MASK = rh.mask,
  IMAGES_LIST = 'rh.T1.txt',
  IMAGING_NAME = 'rh.T1.nii.gz',
  to_gifti = rh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#left
#results.T1.reliability.lh =
doit(
  DATADIR,
  reliability_depth,
  'tests/T1/reliability.lh',
  MASK = lh.cortex.mask,
  IMAGES_LIST = 'lh.T1.txt',
  IMAGING_NAME = 'lh.T1.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

#left
#results.T1.reliability.lh =
doit(
  DATADIR,
  reliability_depth,
  'tests/T1/reliability-mask.lh',
  MASK = lh.mask,
  IMAGES_LIST = 'lh.T1.txt',
  IMAGING_NAME = 'lh.T1.nii.gz',
  to_gifti = lh.gii,
  alpha = alphavertex,
  NPERMS = NPERMS,
  shuffle_by = shuffle_by
)

}

