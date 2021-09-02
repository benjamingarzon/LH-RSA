rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

# remove motion outliers...?

# difference alpha  trained > untrained increases with time and then goes to 0. More for intervention group
# check that not explained by noise!!
# use aroma to denoise??
# why is clf accuracy so high everywhere??
# correct for configuration x time
# correct for FD. FD vs type of 
# correct for number of correct trials
# do it for tessellation
WIDTH = 30
HEIGHT = 24
DPI = 1000
setwd("~/Software/LeftHand/analysis")
source('./plot_funcs.R')
source('./load_covariates.R')

meas = 'correlation'
#meas = 'cosine'
#meas = 'xnobis'
suffix0 = 'mask-cross-perm'
suffix1 = '-same'
suffix = paste0(suffix0, suffix1)
#suffix = 'mask-cross'
sem = function(x) sd(x, na.rm = T)/sqrt(length(x))
data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_', suffix0, '.csv')
#data_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_mask-cross.csv'
#output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/clf_acc.csv'
output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/alpha_diff.csv'
figs_dir = '/data/lv0/MotorSkill/figs/variability'
data = read.table(data_file, header = T, sep = ',') %>% arrange(subject, session, hemi, label)

data$GROUP = "Intervention"
data$GROUP[grep("sub-lue.2", data$subject)] = "Control"

if (meas == 'alpha_PCM'){
  mymeasures = c("alpha_trained_PCM", "alpha_untrained_PCM")
  data$value = data$alpha_trained_PCM - data$alpha_untrained_PCM
} 

if (meas == 'alpha'){
  mymeasures = c("alpha_trained", "alpha_untrained")
  data$value = data$alpha_trained - data$alpha_untrained
} 

if (meas == 'clf'){
  mymeasures = c("clf_acc_trained", "clf_acc_untrained", "clf_acc_trained_untrained")
  data$value = data$clf_acc_trained - data$clf_acc_untrained
  ylabel = "Classification accuracy"
  ylimit = c(.4, .7)
  ylimit_diff = c(-.1, .1)
} 

if (meas == 'euclidean'){
  mymeasures = c("euclidean_trained_different", "euclidean_untrained_different", "euclidean_trained_untrained")
  data$value = data$euclidean_trained_different - data$euclidean_untrained_different 
}

if (meas == 'xnobis'){
  mymeasures = c("xnobis_trained_same", "xnobis_untrained_same", "xnobis_trained_different", "xnobis_untrained_different", "xnobis_trained_untrained")
  mymeasures = c("xnobis_trained_different", "xnobis_untrained_different", "xnobis_trained_untrained")
  data$value = data$xnobis_untrained_same - data$xnobis_trained_same 
  ylabel = "Cross-nobis distance"
  ylimit = c(-.01, .07)
  ylimit_diff = c(-.01, .01)
} 

if (meas == 'cosine'){
  mymeasures = c("cosine_trained_same", "cosine_untrained_same")
#  mymeasures = c("cosine_trained_same", "cosine_untrained_same", "cosine_trained_different", "cosine_untrained_different", "cosine_trained_untrained")
  data$value = data$cosine_trained_same - data$cosine_untrained_same 
  ylabel = "Cosine distance"
  ylimit = c(.8, 1)
  ylimit_diff = c(-.01, .02)
  
  } 

if (meas == 'correlation'){
#   mymeasures = c("correlation_same", "correlation_different")
#  mymeasures = c("correlation_trained_same", "correlation_untrained_same", "correlation_trained_different", 
   #                 "correlation_untrained_different", "correlation_trained_untrained")
  mymeasures = c("correlation_trained_same", "correlation_untrained_same")
#  mymeasures = c("correlation_trained_different", "correlation_untrained_different")
  data$value = data$correlation_trained_same - data$correlation_untrained_same 

#  data$value = data$correlation_trained_different - data$correlation_untrained_different 
  ylabel = "Correlation distance"
  ylimit = c(.8, 1)
  ylimit_diff = c(-.005, .015)
} 

if (meas == 'mean_signal'){
  mymeasures = c("mean_signal_trained", "mean_signal_untrained")

  data$value = data$mean_signal_trained - data$mean_signal_untrained 
} 

##mymeasures = c("clf_acc_trained", "clf_acc_untrained")
#c("mean_signal_trained", "mean_signal_untrained") #c("alpha_trained", "alpha_untrained")
#mymeasures = c("crossnobis_trained", "crossnobis_untrained") #c("mean_signal_trained", "mean_signal_untrained") #c("alpha_trained", "alpha_untrained")


# save it 
output.data = data%>%
  group_by(label)%>%
  summarise(value=mean(value)) %>% arrange(-value) #, value_perm = mean(value_perm)

write.table(output.data, file = output_file, sep = ',', col.names = T, row.names = F)
#stophere

nrois = 5
# check data
plot(sort(table(data$subject))/nrois/2, las = 2) # how many sessions per subject
plot(sort(table(data$session))/nrois/2, las = 2) # how many sessions per timepoint

# relabel rois and separate hemispheres
control_labels = c("R_C1", "L_C1")
#data = data %>% filter( !label %in% control_labels)
incomplete_subjects = c("sub-lue5207", "sub-lue3202", "sub-lue1201")
data$label  = gsub("R_", "Right ", data$label)
data$label  = gsub("L_", "Left ", data$label)
data$label  = gsub("C1", "Control Region", data$label)
data$uni_label  = gsub("Right ", "Average ", data$label)
data$uni_label  = gsub("Left ", "Average ", data$uni_label)

data = data %>% filter(!is.infinite(value))  %>% filter(!subject %in% incomplete_subjects, valid_runs >3) 

data = data %>% group_by(subject) %>% 
  mutate(value=markoutliersIQR(value)) %>% 
  filter(!is.na(value)) %>% ungroup()

data = merge(data, covars.table %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
             by = c("subject", "session", "GROUP"))

data = left_join(data, motion %>% group_by(SUBJECT, TP) %>% 
                   summarise(FD = mean (FD)) %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
                 by = c("subject", "session"))

data.melt = melt(data, 
                 id.vars = c("subject", "session", "GROUP", "hemi", "label", 
                             "FD", "SYSTEM", "CONFIGURATION", "TRAINING", "TRAINING.Q"),
                 variable.name = "MEASURE", 
                 value.name = "value") %>% mutate(SUBJECT = subject, value = as.numeric(value)) %>%
  filter(MEASURE %in% mymeasures) 

# relabel measures 


# model the data + as.factor(CONFIGURATION)*as.factor(session)
#

if (F) {
data.melt$value.corr  = 0

# for (l in unique(data.melt$label)) {
#   print(l)
#   model = lmer(clf_acc ~ 1 + FD + SYSTEM + (1 |subject), data = data %>% filter(label == l))
#   print(summary(model))
# }
for (l in unique(data.melt$label)) {
  print(l)
model = lmer(value ~ 1 + FD + SYSTEM + as.factor(CONFIGURATION) + GROUP*MEASURE*(TRAINING + TRAINING.Q) + 
                 + (1 + TRAINING + TRAINING.Q|SUBJECT), data = data.melt %>% filter(label == l))
print(summary(model))
data.melt$value.corr[data.melt$label == l] = predict(model, data.melt %>% filter(label == l) %>% 
                                                           mutate(FD = 0,
                                                           SYSTEM = 'Classic',
                                                           CONFIGURATION = '1'
                                                           )) + resid(model)
}

} else {
  data.melt$value.corr  = data.melt$value
  
}
# myplot = ggplot(data, aes(
#   x = session,
#   y = value,
#   col = group,
#   group = subject
# )) + geom_line(size = 0.2, alpha = 0.5) + 
#   facet_wrap(. ~ label, ncol = 5) + 
#   theme_lh() + scale_colour_manual(values = myPalette)

#print(myplot)

# average

data.mean = data%>%
  group_by(session, GROUP, label)%>%summarise(val.mean = mean(value, na.rm = T), 
                                       val.sem = sem(value)) 
data.uni = data %>%
  group_by(subject, session, GROUP, uni_label)%>%summarise(value = mean(value, na.rm = T)) %>%
  group_by(session, GROUP, uni_label)%>%summarise(val.mean = mean(value, na.rm = T), 
                                                         val.sem = sem(value)) %>% mutate(label = uni_label)
data.mean = rbind(data.mean, data.uni)

myplot = ggplot(data.mean , aes( #%>% filter(hemi == 'rh')
  x = session,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = GROUP,
  group = GROUP
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_wrap(. ~ label, ncol = nrois) + 
  ylab(ylabel) +
  xlab('Test session') + 
  ylim(ylimit_diff) + 
  theme_lh() + 
  scale_colour_manual(values = myPalette) +
  theme(legend.title = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #legend.position = 'bottom', 

print(myplot)


# hemi
clean_measure_names = function(x){
  
    sapply(x, function(z) {
    z = gsub('acc_', '', z)
    str_to_title(paste(strsplit(as.character(z), '_')[[1]][-1], collapse = ' '))})
}
data.melt = data.melt%>% mutate(MEASURE = clean_measure_names(MEASURE))
data.mean = data.melt%>%
  group_by(session, GROUP, label, MEASURE)%>%summarise(val.mean = mean(value.corr, na.rm = T), 
                                                    val.sem = sem(value.corr))%>% 
  mutate(group_measure = paste(GROUP, MEASURE)) 

myplot.all = ggplot(data.mean, aes(
  x = session,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = MEASURE,
  group = group_measure,
  linetype = GROUP
  )) + 
  geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_wrap(. ~label, ncol = nrois) + 
  ylab(ylabel) +
  xlab('Test session') + 
  ylim(ylimit) + 
  theme_lh() + 
  scale_colour_manual(values = myPalette) + 
  theme(text = element_text(size = 8), #legend.position = 'bottom', 
        legend.title = element_blank(),
        legend.box = 'vertical',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(myplot.all)

ggsave(
  file.path(figs_dir, paste(meas, suffix, 'difference.png', sep = '-')),
  plot = myplot,
  width = WIDTH,
  height = HEIGHT,
  dpi = DPI,
  units = 'cm'
)
ggsave(
  file.path(figs_dir, paste(meas, suffix, 'all.png', sep = '-')),
  plot = myplot.all,
  width = WIDTH,
  height = HEIGHT,
  dpi = DPI,
  units = 'cm'
)

#alpha
#plot(data$alpha_trained, data$alpha_untrained, pch = 20, cex = 0.2)
abline(0, 1)
alphas = dplyr::select(data, c(alpha_trained, alpha_untrained, alpha_trained_PCM, alpha_untrained_PCM))
alphas = alphas[rowMeans(alphas) > -Inf, ]
alphas = alphas[alphas$alpha_untrained_PCM > -2, ]
cor(alphas[complete.cases(alphas), ])
cor.test(data$alpha_trained, data$alpha_untrained)
cor.test(data$alpha_trained[!is.infinite(data$alpha_trained_PCM)], data$alpha_trained_PCM[!is.infinite(data$alpha_trained_PCM)])
cor.test(data$alpha_untrained[!is.infinite(data$alpha_untrained_PCM)], data$alpha_untrained_PCM[!is.infinite(data$alpha_untrained_PCM)])
cor.test(data$alpha_trained[!is.infinite(data$alpha_trained_PCM)], data$alpha_trained_PCM[!is.infinite(data$alpha_trained_PCM)])

#plot(data$alpha_untrained, data$alpha_untrained_PCM, ylim = c(-2, 4))

# myplot = ggplot(data%>% filter( !label %in% control_labels) %>% group_by(session, group, hemi, CONFIGURATION) %>% 
#                   summarise(val.mean = mean(value), val.sem = sem(value)), aes(
#   x = session,
#   y = val.mean,
#   ymin = val.mean - val.sem, 
#   ymax = val.mean + val.sem, 
#   col = group,
#   group = group
# )) + geom_line() + 
#   geom_point() + 
#   geom_errorbar() + 
#   facet_grid(CONFIGURATION ~ hemi) + 
#   #ylim(0, 1) + 
#   theme_lh() + scale_colour_manual(values = myPalette)
# 
# print(myplot)
