rm(list = ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(pracma)
library(freesurfer)
library(caret)
library(glmnet)
library(forecast)
library(e1071)
library(doParallel)
library(foreach)


# correct for system used
#set.seed(1)
#sync permutations!!!!!!!!!

# can put back system as covariate
# project to sphere
# pca
# fit models
# null distribution

# make sure variances are ok xx
# project weights xx
# iterate correlations
# make sure that null distribution is correct
# remove 7 timepoint
# screeplot
# RBM
# remove outliers

# try on whole brain xx
# try repeated cv xx
# output results xx
# balance classes xx
# cross-validation... xx
#tuning xx
# permutation xx
# remove effect of classic... xx
# parallelize xx

cat("\014")

source('~/Software/ImageLMMR/ImageLMMR.R')
setwd("~/Software/LeftHand/analysis")
source('./structural_analysis_funcs.R')

SUBJECTS_DIR = '~/Data/LeftHand/Lund1/freesurfer/'

OUTPUTDIR = 'multivariate-tess'
NPROCS = 38
NINDS = NPROCS * 5 # how often to dump the data

source('./load_covariates.R')

###############################################################
# Define the pipeline
###############################################################
doit = structural_analysis_multivariate

###############################################################
# Do it
###############################################################
vol.mask = 'mask.nii.gz'
rh.cortex.mask = 'rh.cortex.mask.nii.gz'
lh.cortex.mask = 'lh.cortex.mask.nii.gz'

shuffle_by = c('BETWEEN', 'WITHIN')

alphavoxel = 0.05
alphavertex = 0.025

NPERMS = 200 

#####################
# VBM
#####################

#DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/cat12crossbias8_10'
DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/cat12'

annotation.MNI = c('/home/benjamin.garzon/Data/LeftHand/Lund1/labels/cvs_avg35_inMNI152/Icosahedron-162/Icosahedron-162.aparc.cat.nii.gz',
                   '/home/benjamin.garzon/Data/LeftHand/Lund1/labels/cvs_avg35_inMNI152/Icosahedron-162/LUT.txt')
                   
labels.MNI = '/home/benjamin.garzon/Software/LeftHand/masks/MNI.tessellation162.txt'

if (F){
results.multivariate.vol = doit(
  DATADIR,
  MASK = vol.mask,
  IMAGES_LIST = 'image_list.txt',
  IMAGING_NAME = 'data.nii.gz',
  HEMI = NULL,
  ANNOT_FILE = annotation.MNI,
  LABEL_FILE = labels.MNI,
  ACC_FILE = 'MNI.accuracy',
  NPERMS = NPERMS,
  OUTDIR = file.path(DATADIR, 'tests', OUTPUTDIR)
) 

save(
  results.multivariate.vol,
  file = file.path(DATADIR, 'tests', OUTPUTDIR, 'results.rds')
)
}

#####################
# thickness
#####################

DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results_unsmoothed'
#annotation.rh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations/GlasserParc/rh.HCP-MMP1.annot'
#labels.rh = '/home/benjamin.garzon/Software/LeftHand/masks/rh.parcels.txt'
#annotation.lh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations/GlasserParc/lh.HCP-MMP1.annot'
#labels.lh = '/home/benjamin.garzon/Software/LeftHand/masks/lh.parcels.txt'

annotation.rh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations/fs_LR_32/Icosahedron-162.fsaverage.R.annot'
labels.rh = '/home/benjamin.garzon/Software/LeftHand/masks/rh.tessellation162.txt'
annotation.lh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations//fs_LR_32/Icosahedron-162.fsaverage.L.annot'
labels.lh = '/home/benjamin.garzon/Software/LeftHand/masks/lh.tessellation162.txt'

results.multivariate.rh = doit(
  DATADIR,
  MASK = rh.cortex.mask,
  IMAGES_LIST = 'rh.thickness.txt',
  IMAGING_NAME = 'rh.thickness.nii.gz',
  HEMI = 'rh',
  ANNOT_FILE = annotation.rh,
  LABEL_FILE = labels.rh,
  ACC_FILE = 'rh.accuracy',
  NPERMS = NPERMS,
  OUTDIR = file.path(DATADIR, 'tests', OUTPUTDIR)
) # remove outliers

results.multivariate.lh = doit(
  DATADIR,
  MASK = lh.cortex.mask,
  IMAGES_LIST = 'lh.thickness.txt',
  IMAGING_NAME = 'lh.thickness.nii.gz',
  HEMI = 'lh',
  ANNOT_FILE = annotation.lh,
  LABEL_FILE = labels.lh,
  ACC_FILE = 'lh.accuracy',
  NPERMS = NPERMS,
  OUTDIR = file.path(DATADIR, 'tests', OUTPUTDIR)
) # remove outliers

save(
  results.multivariate.lh,
  results.multivariate.rh,
  file = file.path(DATADIR, 'tests', OUTPUTDIR, 'results.rds')
)
