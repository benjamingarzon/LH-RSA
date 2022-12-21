# compose figs for structural analyses

rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggpubr)
library(png)
library(lme4)
library(lmerTest)
library(reshape2)
library(reshape)
library(doParallel)
library(foreach)

setwd("~/Software/LeftHand/analysis")

FIGS_DIR = '~/Data/LeftHand/Lund1/figs/structure'
DPI = 1000
FONT.SIZE = 14
NPROCS = 35
NSAMPLES = 500
NVOXELS = 200
N = 30 # subjects per group

maxvals = seq(0, 5, 1) / 100

anisotropic = F

# Vertex
isROI = 0
alpha = 0.001

mytitle = 'Grey matter volume, voxel-wise analysis'
doGMV = T
#source("./do_change_simulations.R")

mytitle = 'Cortical thickness, vertex-wise analysis'
doGMV = F
#source("./do_change_simulations.R")

#ROI
alpha = 0.05
doGMV = T

mytitle = 'Grey matter volume, ROI analysis, 33%'
isROI = .33
source("./do_change_simulations.R")

mytitle = 'Grey matter volume, ROI analysis, 66%'
isROI = .66
source("./do_change_simulations.R")

mytitle = 'Grey matter volume, ROI analysis, 100%'
isROI = 1
source("./do_change_simulations.R")

doGMV = F

mytitle = 'Cortical thickness, ROI analysis, 33%'
isROI = .33
source("./do_change_simulations.R")

mytitle = 'Cortical thickness, ROI analysis, 66%'
isROI = .66
source("./do_change_simulations.R")

mytitle = 'Cortical thickness, ROI analysis, 100%'
isROI = 1
source("./do_change_simulations.R")

#img.thickness.anisotropic = readPNG(file.path(FIGS_DIR, 'PowerThicknessAnisotropic.png'))
#img.thickness.anisotropic = ggplot() + background_image(img.thickness.anisotropic)

img.GMV = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerGMV.png')))
img.thickness.isotropic = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerThicknessIsotropic.png')))

#img.GMV.ROI = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerGMV-ROI.png')))
#img.thickness.isotropic.ROI = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerThicknessIsotropic-ROI.png')))

img.GMV.ROI.1 = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerGMV-ROI-0.33.png')))
img.thickness.isotropic.ROI.1 = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerThicknessIsotropic-ROI-0.33.png')))

img.GMV.ROI.2 = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerGMV-ROI-0.66.png')))
img.thickness.isotropic.ROI.2 = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerThicknessIsotropic-ROI-0.66.png')))

img.GMV.ROI.3 = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerGMV-ROI-1.png')))
img.thickness.isotropic.ROI.3 = ggplot() + background_image(readPNG(file.path(FIGS_DIR, 'PowerThicknessIsotropic-ROI-1.png')))

plot.final = ggarrange(
  img.thickness.isotropic, 
  img.GMV, 
  img.thickness.isotropic.ROI.1, 
  img.GMV.ROI.1, 
  img.thickness.isotropic.ROI.2, 
  img.GMV.ROI.2, 
  img.thickness.isotropic.ROI.3, 
  img.GMV.ROI.3, 
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"), 
  ncol = 2, nrow = 4, font.label = list(size = FONT.SIZE))

ggsave(file.path(FIGS_DIR, 'Power.png'), plot = plot.final, dpi = DPI, width = 16, height = 26, units = 'cm')
  

