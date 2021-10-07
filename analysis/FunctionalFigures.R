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
HEIGHT.1 = 10
WIDTH = 14
HEIGHT.2 = 20
HEIGHT.3 = 30
FONT.LABEL.W = list(size = 14, color = "white", face = "bold", family = NULL)
FONT.LABEL.B = list(size = 14, color = "black", face = "bold", family = NULL)

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
                      labels = c("A", "B", "C"), ncol = 3, nrow = 1, font.label = FONT.LABEL.W)

act.plot = ggarrange(act.plot1, activation_D, labels = c("", "D"), font.label = FONT.LABEL.B, 
                      ncol = 1, nrow = 2)

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
                         ncol = 1, nrow = 2)

ggsave(filename = file.path(FIG_DIR, "Fig_UntrainedSameDiffSupp.png"), 
       plot = samediff.plotsupp,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.1)


#distances_imgA = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-comparison_prereg-Best-map.png'))
distances_imgA = readPNG(file.path(FIG_DIR, 'variability/xnobis-mask-cross-runprew-same-all.png'))
distances_imgB = readPNG(file.path(FIG_DIR, 'variability/xnobis-mask-cross-runprew-different-all.png'))
distances_imgC = readPNG(file.path(FIG_DIR, 'variability/xcorrelation-mask-cross-runprew-same-all.png'))
distances_imgD = readPNG(file.path(FIG_DIR, 'variability/xcorrelation-mask-cross-runprew-different-all.png'))

distances_A <- ggplot() + background_image(distances_imgA) 
distances_B <- ggplot() + background_image(distances_imgB) 
distances_C <- ggplot() + background_image(distances_imgC) 
distances_D <- ggplot() + background_image(distances_imgD) 

#distances_A = ggarrange(distances_A,
#                          labels = c("A"), ncol = 1, nrow = 1, font.label = FONT.LABEL.W)

distances.plot = ggarrange(distances_A, distances_B, distances_C, distances_D,
                          labels = c("", "B", "C", "D"), ncol = 2, nrow = 2)

ggsave(filename = file.path(FIG_DIR, "Fig_PatternDistances.png"), 
       plot = distances.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT.2)



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