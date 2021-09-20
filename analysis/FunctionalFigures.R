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
HEIGHT = 12
WIDTH = 10
FONT.LABEL.W = list(size = 14, color = "white", face = "bold", family = NULL)
FONT.LABEL.B = list(size = 14, color = "black", face = "bold", family = NULL)

########################
# Activation Figures
########################
activation_imgA = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_whole_prereg-INTERCEPT_p_fdr-map.png'))
activation_imgB = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_prereg_groupxcondition-Omni_p_fdr-map.png'))
activation_imgC = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_prereg_groupxconditionxtraining-Omni_p_fdr-map.png'))
activation_imgD = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_prereg_groupxcondition-Omni_p_fdr-data.png'))
activation_imgD = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_prereg_groupxconditionxtraining-Omni_p_fdr-data.png'))

activation_A <- ggplot() + background_image(activation_imgA) 
activation_B <- ggplot() + background_image(activation_imgB) 
activation_C <- ggplot() + background_image(activation_imgC) 
activation_D <- ggplot() + background_image(activation_imgD) 
activation_E <- ggplot() + background_image(activation_imgE) 

########################
# Session effects
########################
session_imgA = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_whole_prereg-INTERCEPT_p_fdr-map.png'))
session_imgB = readPNG(file.path(FIG_DIR, 'Trained_Untrained/tests-quadratic_prereg_groupxcondition-Omni_p_fdr-map.png'))

session_A <- ggplot() + background_image(session_imgA) 
session_B <- ggplot() + background_image(session_imgB) 




########################
# Combine figures
########################

# activation
act.plot1 = ggarrange(activation_A, activation_B, activation_C,
                      labels = c("A", "B", "C"), ncol = 3, nrow = 1, font.label = FONT.LABEL.W)
act.plot2 = ggarrange(activation_D,
                      labels = "D", ncol = 1, nrow = 1, font.label = FONT.LABEL.B)
act.plot2 = ggarrange(activation_E,
                      labels = "E", ncol = 1, nrow = 1, font.label = FONT.LABEL.B)

act.plot = ggarrange(act.plot1, act.plot2, act.plot3,
                      ncol = 1, nrow = 3, heights = c(1, 3))

ggsave(filename = file.path(FIG_DIR, "Fig_ActivationEffects.png"), 
       plot = act.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT)


# session effects
sess.plot = ggarrange(sess.plot1, sess.plot2,
                     ncol = 1, nrow = 2, heights = c(1, 2))

ggsave(filename = file.path(FIG_DIR, "Fig_ActivationSessionEffects.png"), 
       plot = act.plot,
       dpi = DPI, 
       width = WIDTH, 
       height = HEIGHT)


# multivariate patterns