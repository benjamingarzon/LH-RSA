# compose figs for structural analyses

rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggpubr)
library(png)

FIGS_DIR='~/Data/LeftHand/Lund1/figs/structure'
DPI = 1000
FONT.SIZE = 14

######################################
# Reliability
######################################

#####################
# Plot reliability histograms
#####################
DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'
load(file.path(DATADIR, 'tests', 'thickness', 'reliability-mask.rh', 'results.rda'))
stats = results$stats
load(file.path(DATADIR, 'tests', 'thickness', 'reliability-mask.lh', 'results.rda'))
stats.thickness = c(stats, results$stats)
print(paste(median(stats.thickness), sd(stats.thickness)))
reliability.thickness.plot = qplot(stats.thickness, geom="histogram", binwidth = .02, xlab= 'ICC', 
                                   main = "Cortical thickness", xlim = c(0, 1)) + theme_classic() 


# T1
load(file.path(DATADIR, 'tests', 'T1', 'reliability-mask.rh', 'results.rda'))
stats = results$stats
load(file.path(DATADIR, 'tests', 'T1', 'reliability-mask.lh', 'results.rda'))
stats.T1 = c(stats, results$stats)
print(paste(median(stats.T1), sd(stats.T1)))
reliability.T1.plot = qplot(stats.T1, geom="histogram", binwidth = .02, xlab= 'ICC', 
                                   main = "T1 relaxation time", xlim = c(0, 1)) + theme_classic() 

DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/cat12'
load(file.path(DATADIR, 'tests', 'reliability-mask', 'results.rda'))
print(paste(median(results$stats), sd(results$stats)))
reliability.VBM.plot = qplot(results$stats, geom="histogram", binwidth = .02, xlab= 'ICC', 
                             main = "Grey matter volume", xlim = c(0,1)) + theme_classic() 

img.thickness = readPNG(file.path(FIGS_DIR, 'reliability-thickness.png'))
img.VBM = readPNG(file.path(FIGS_DIR, 'reliability-VBM-FSL.png'))
img.T1 = readPNG(file.path(FIGS_DIR, 'reliability-T1.png'))

img.thickness = ggplot() + background_image(img.thickness)
img.VBM = ggplot() + background_image(img.VBM)
img.T1 = ggplot() + background_image(img.T1)

img.thickness = ggarrange(img.thickness, labels = "A", font.label = list(size = FONT.SIZE, color = "white"))
img.VBM = ggarrange(img.VBM, labels = "B", font.label = list(size = FONT.SIZE, color = "white"))
img.T1 = ggarrange(img.T1, labels = "C", font.label = list(size = FONT.SIZE, color = "white"))


reliability.plot = ggarrange(img.thickness, reliability.thickness.plot, 
                             img.VBM, reliability.VBM.plot, 
                             img.T1, reliability.T1.plot, 
                             labels = c("", "D", "", "E", "", "F"), 
                          ncol = 2, nrow = 3, font.label = list(size = FONT.SIZE))

ggsave(file.path(FIGS_DIR, 'Reliability-analyses.png'), plot = reliability.plot, dpi = DPI, width = 11, height = 9)

stophere
######################################
# Univariate analyses
######################################

######################################
# Multivariate analyses
######################################

img.map.1 = readPNG(file.path(FIGS_DIR, 'multivar-thickness', 'multivar-tess-thickness.png'))
img.accu.1 = readPNG(file.path(FIGS_DIR, 'multivar-thickness', 'Accuracy.png'))
img.regions.1 = readPNG(file.path(FIGS_DIR, 'multivar-thickness', 'IndividualRegions.png'))
img.timecourse.1 = readPNG(file.path(FIGS_DIR, 'multivar-thickness', 'Timecourse.png'))

img.map.2 = readPNG(file.path(FIGS_DIR, 'multivar-VBM', 'multivar-tess-vbm.png'))
img.accu.2 = readPNG(file.path(FIGS_DIR, 'multivar-VBM', 'Accuracy.png'))
img.regions.2 = readPNG(file.path(FIGS_DIR, 'multivar-VBM', 'IndividualRegions.png'))
img.timecourse.2 = readPNG(file.path(FIGS_DIR, 'multivar-thickness', 'Timecourse.png')) ## change

img.map.1 = ggplot() + background_image(img.map.1)
img.accu.1 = ggplot() + background_image(img.accu.1)
img.regions.1 = ggplot() + background_image(img.regions.1)
img.timecourse.1 = ggplot() + background_image(img.timecourse.1)

img.map.2 = ggplot() + background_image(img.map.2)
img.accu.2 = ggplot() + background_image(img.accu.2)
img.regions.2 = ggplot() + background_image(img.regions.2)
img.timecourse.2 = ggplot() + background_image(img.timecourse.2)

multivar.plot.thickness = ggarrange(img.map.1, img.regions.1, img.accu.1, img.timecourse.1,
                     labels = c("A", "B", "C", "D"), 
                     ncol = 2, nrow = 2, font.label = list(size = FONT.SIZE))


ggsave(file.path(FIGS_DIR, 'Multivariate-thickness-analyses.png'), dpi = DPI, width = 6.5, height = 4)

multivar.plot.VBM = ggarrange(img.map.2, img.accu.2, img.regions.2, img.timecourse.2,
                                    labels = c("A", "C", "B", "D"), 
                                    ncol = 2, nrow = 2, font.label = list(size = FONT.SIZE))

ggsave(file.path(FIGS_DIR, 'Multivariate-VBM-analyses.png'), dpi = DPI, width = 6, height = 3)
