# Create figure explaining design
rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggpubr)
library(png)

setwd("~/Software/LeftHand/analysis")
source("./plot_funcs.R")
FIG_DIR = '~/Data/LeftHand/Lund1/figs/'

DPI = 1000
WIDTH = 16
HEIGHT.1 = 8 # 1 row
HEIGHT.1a = 5 # 1 row
HEIGHT.2 = 18 # 2 rows
HEIGHT.3 = 22 # 3 rows
FONT.LABEL.W = list(size = 14, color = "white", face = "bold", family = NULL)
FONT.LABEL.W2 = list(size = 20, color = "white", face = "bold", family = NULL)
FONT.LABEL.B = list(size = 14, color = "black", face = "bold", family = NULL)
FONT.LABEL.B2 = list(size = 20, color = "black", face = "bold", family = NULL)

########################
# Activation Figures
########################
activation_imgA = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-asymptotic_whole_prereg-INTERCEPT_p_fdr-map.png'))
activation_imgB = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-comparison_prereg-Best-map.png'))
activation_imgC = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-asymptotic_prereg_groupxconditionxtraining-Omni_p_fdr-map.png'))
activation_imgD = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-asymptotic_prereg_groupxconditionxtraining-Omni_p_fdr-data.png'))
activation_imgE = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-asymptotic_prereg_groupxtraining-Omni_p_fdr-map.png'))
activation_imgF = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-asymptotic_prereg_groupxtraining-Omni_p_fdr-data.png'))

activation_A <- ggplot() + background_image(activation_imgA) 
activation_B <- ggplot() + background_image(activation_imgB) 
activation_C <- ggplot() + background_image(activation_imgC) 
activation_D <- ggplot() + background_image(activation_imgD) 
activation_E <- ggplot() + background_image(activation_imgE) 
activation_F <- ggplot() + background_image(activation_imgF) 

# activation
act.plot1 = ggarrange(activation_A, activation_B, activation_C,
                      labels = c("A", "B", "C"), ncol = 3, nrow = 1, font.label = FONT.LABEL.W2)

act.plot = ggarrange(act.plot1, activation_D, labels = c("", "D"), font.label = FONT.LABEL.B, 
                      ncol = 1, nrow = 2, heights = c(1, 0.9))

ggsave(filename = file.path(FIG_DIR, "Fig_ActivationEffects.png"), 
       plot = act.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.1)

act.plotsupp1 = ggarrange(activation_E,
                      labels = c("A"), ncol = 1, nrow = 1, font.label = FONT.LABEL.W)
act.plotsupp2 = ggarrange(activation_F,
                          labels = c("B"), ncol = 1, nrow = 1, widths = 1, font.label = FONT.LABEL.B)
act.plotsupp = ggarrange(act.plotsupp1, act.plotsupp2,
                          ncol = 2, nrow = 1, widths = c(1, 4))

ggsave(filename = file.path(FIG_DIR, "Fig_ActivationEffectsSupp.png"), 
       plot = act.plotsupp,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.1)


########################
# Multivariate patterns
########################
# variability 
alpha_img = readPNG(file.path(FIG_DIR, 'variability/alpha-mask-cross-runprew-perm-all-all.png'))

alpha.plot <- ggplot() + background_image(alpha_img) 


ggsave(filename = file.path(FIG_DIR, "Fig_Alpha.png"), 
       plot = alpha.plot,
       dpi = DPI) 
       #width = WIDTH, 
       #height = HEIGHT.1)



# same vs diff
samediff_imgA = readPNG(file.path(FIG_DIR, 'variability/xnobis-mask-cross-runprew-untrained-difference.png'))
samediff_imgB = readPNG(file.path(FIG_DIR, 'variability/xcorrelation-mask-cross-runprew-untrained-difference.png'))
samediff_imgC = readPNG(file.path(FIG_DIR, 'variability/xnobis-mask-cross-runprew-perm-untrained-difference.png'))
samediff_imgD = readPNG(file.path(FIG_DIR, 'variability/xcorrelation-mask-cross-runprew-perm-untrained-difference.png'))

samediff_A <- ggplot() + background_image(samediff_imgA) 
samediff_B <- ggplot() + background_image(samediff_imgB) 
samediff_C <- ggplot() + background_image(samediff_imgC) 
samediff_D <- ggplot() + background_image(samediff_imgD) 

samediff.plotsupp = ggarrange(samediff_A, samediff_B,
                              labels = c("A", "B"), 
                         ncol = 1, nrow = 2, font.label = FONT.LABEL.B2)

ggsave(filename = file.path(FIG_DIR, "Fig_UntrainedSameDiffSupp.png"), 
       plot = samediff.plotsupp,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.2)

# distances
distances_imgA = readPNG(file.path(FIG_DIR, 'MultivarROIs.png'))
distances_imgB = readPNG(file.path(FIG_DIR, 'variability/xnobis-mask-cross-runprew-same-all.png'))
distances_imgC = readPNG(file.path(FIG_DIR, 'variability/xnobis-mask-cross-runprew-different-all.png'))
distances_imgD = readPNG(file.path(FIG_DIR, 'variability/xcorrelation-mask-cross-runprew-same-all.png'))
distances_imgE = readPNG(file.path(FIG_DIR, 'variability/xcorrelation-mask-cross-runprew-different-all.png'))

distances_A <- ggplot() + background_image(distances_imgA) 
distances_B <- ggplot() + background_image(distances_imgB) 
distances_C <- ggplot() + background_image(distances_imgC) 
distances_D <- ggplot() + background_image(distances_imgD) 
distances_E <- ggplot() + background_image(distances_imgE) 

distances_1 = ggarrange(distances_A,
                          labels = c("A"), ncol = 1, nrow = 1, font.label = FONT.LABEL.W2)
distances_2 = ggarrange(distances_B, distances_C,
                        labels = c("B", "C"), ncol = 2, nrow = 1, font.label = FONT.LABEL.B2)
distances_3 = ggarrange(distances_D, distances_E,
                        labels = c("D", "E"), ncol = 2, nrow = 1, font.label = FONT.LABEL.B2)

distances.plot = ggarrange(distances_1, distances_2, distances_3,
                  ncol = 1, nrow = 3, heights = c(1, 1.8, 1.8))

ggsave(filename = file.path(FIG_DIR, "Fig_PatternDistances.png"), 
       plot = distances.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.3)

# classification
clf_imgA = readPNG(file.path(FIG_DIR, 'variability/clf-mask-cross-runprew-all-all.png'))
clf_imgB = readPNG(file.path(FIG_DIR, 'variability/clf-mask-cross-runprew-perm-all-all.png'))

clf_imgA <- ggplot() + background_image(clf_imgA) 
clf_imgB <- ggplot() + background_image(clf_imgB) 

clf.plot = ggarrange(clf_imgA, clf_imgB, labels = c("A", "B"), ncol = 2, nrow = 1)

ggsave(filename = file.path(FIG_DIR, "Fig_ClassifierAccuracy.png"), 
       plot = clf.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.1)

clf_imgA = readPNG(file.path(FIG_DIR, 'variability/clf-mask-cross-runprew-different-all.png'))
clf_imgB = readPNG(file.path(FIG_DIR, 'variability/clf-mask-cross-runprew-perm-different-all.png'))

clf_imgA <- ggplot() + background_image(clf_imgA) 
clf_imgB <- ggplot() + background_image(clf_imgB) 

clf.plot = ggarrange(clf_imgA, clf_imgB, labels = c("A", "B"), ncol = 2, nrow = 1)

ggsave(filename = file.path(FIG_DIR, "Fig_ClassifierAccuracyDifferent.png"), 
       plot = clf.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.1)

# tessellation
tess_imgA = readPNG(file.path(FIG_DIR, 'tessellation/GROUPxTRAINING_p_uncorrected.png'))
tess_imgB = readPNG(file.path(FIG_DIR, 'tessellation/GROUPxTRAINING_p.png'))

tess_imgA <- ggplot() + background_image(tess_imgA) 
tess_imgB <- ggplot() + background_image(tess_imgB) 

tess.plot = ggarrange(tess_imgA, tess_imgB, labels = c("A", "B"), ncol = 2, nrow = 1, 
                      font.label = FONT.LABEL.W2)

ggsave(filename = file.path(FIG_DIR, "Fig_Tessellation.png"), 
       plot = tess.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.1a)


stophere

#WD=/data/lv0/MotorSkill/fmriprep/analysis/surf/
#OVERLAY_L=$WD/lh.tess-cross-runprewxnobis-GROUPIntervention_TRAINING-GAM.map 
#OVERLAY_R=$WD/rh.tess-cross-runprewxnobis-GROUPIntervention_TRAINING-GAM.map 
#THR_L=.95
#THR_H=1
#FIG=/home/xgarzb@GU.GU.SE/Data/LeftHand/Lund1/figs/tessellation/GROUPxTRAINING_p
#SURF_L=/usr/local/freesurfer/7.1.1-1/subjects/fsaverage/surf/lh.inflated
#SURF_R=/usr/local/freesurfer/7.1.1-1/subjects/fsaverage/surf/rh.inflated
#./show_surface.sh $WD $OVERLAY_L $OVERLAY_R $THR_L $THR_H $FIG $SURF_L $SURF_R

########################
# Session effects
########################
#session_imgA = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_whole_prereg-INTERCEPT_p_fdr-map.png'))
#session_imgB = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-asymptotic_prereg_groupxcondition-Omni_p_fdr-map.png'))

#session_A <- ggplot() + background_image(session_imgA) 
#session_B <- ggplot() + background_image(session_imgB) 

# session effects
#sess.plot = ggarrange(sess.plot1, sess.plot2,
#                     ncol = 1, nrow = 2, heights = c(1, 2))

#ggsave(filename = file.path(FIG_DIR, "Fig_ActivationSessionEffects.png"), 
#       plot = sess.plot,
#       dpi = DPI, 
#       width = WIDTH, 
#       height = HEIGHT)


# multivariate patterns