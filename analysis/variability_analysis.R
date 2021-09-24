library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

# remove motion outliers...?
# correct for FD. FD vs type of 
# correct for number of correct trials


do_variability_analysis = function(meas, suffix0, suffix1) {
suffix = paste0(suffix0, suffix1)
data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_', suffix0, '.csv')
#data_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_mask-cross.csv'
#output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/clf_acc.csv'
output_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/', suffix, '.csv')
figs_dir = '/data/lv0/MotorSkill/figs/variability'
data = read.table(data_file, header = T, sep = ',') %>% arrange(subject, session, hemi, label)

data$GROUP = "Intervention"
data$GROUP[grep("sub-lue.2", data$subject)] = "Control"

# remove first wave
data = data[ -grep('sub-lue1', data$subject), ]

if (meas == 'alpha_PCM'){
  mymeasures = c("alpha_trained_PCM", "alpha_untrained_PCM")
  data$value = data$alpha_trained_PCM - data$alpha_untrained_PCM
} 

if (meas == 'alpha'){
  mymeasures = c("alpha_trained", "alpha_untrained")
  data$value = data$alpha_trained - data$alpha_untrained
  ylabel = "Variability score"
  ylimit = c(1.5, 3)
  ylimit_diff = c(-.1, .1)
} 

if (meas == 'clf'){
  mymeasures = c("clf_acc_trained", "clf_acc_untrained", "clf_acc_trained_untrained")
  data$value = data$clf_acc_trained - data$clf_acc_untrained
  ylabel = "Classification accuracy"
  ylimit = c(.4, .7)
  ylimit_diff = c(-.1, .1)
} 

if (meas == 'euclidean'){
  if (suffix1 == '-same') mymeasures = c("euclidean_trained_same", "euclidean_untrained_same")
  if (suffix1 == '-different') mymeasures = c("euclidean_trained_different", "euclidean_untrained_different")
  if (suffix1 == '-all') mymeasures = c("euclidean_trained_same", "euclidean_untrained_same", "euclidean_trained_different", "euclidean_untrained_different", "euclidean_trained_untrained")
  if (suffix1 == '-same') data$value = data$euclidean_trained_same - data$euclidean_untrained_same 
  else data$value = data$euclidean_trained_different - data$euclidean_untrained_different
  ylabel = "Squared euclidean distance"
  ylimit = c(4, 6)
  ylimit_diff = c(-2, 2)
  
  }

if (meas == 'xnobis'){
  if (suffix1 == '-same') mymeasures = c("xnobis_trained_same", "xnobis_untrained_same")
  if (suffix1 == '-different') mymeasures = c("xnobis_trained_different", "xnobis_untrained_different")
  if (suffix1 == '-all') mymeasures = c("xnobis_trained_same", "xnobis_untrained_same", "xnobis_trained_different", "xnobis_untrained_different", "xnobis_trained_untrained")
  if (suffix1 == '-same') data$value = data$xnobis_trained_same - data$xnobis_untrained_same 
  else data$value = data$xnobis_trained_different - data$xnobis_untrained_different

  ylabel = "Cross-nobis distance"
  ylimit = c(-.01, .06)
  ylimit_diff = c(-.02, .02)
} 

if (meas == 'cosine'){
  if (suffix1 == '-same') mymeasures = c("cosine_trained_same", "cosine_untrained_same")
  if (suffix1 == '-different') mymeasures = c("cosine_trained_different", "cosine_untrained_different")
  if (suffix1 == '-all') mymeasures = c("cosine_trained_same", "cosine_untrained_same", "cosine_trained_different", "cosine_untrained_different", "cosine_trained_untrained")
  if (suffix1 == '-same') data$value = data$cosine_trained_same - data$cosine_untrained_same 
  else data$value = data$cosine_trained_different - data$cosine_untrained_different
  ylabel = "Cosine (z-scored)"
  ylimit = c(0, .2)
  ylimit_diff = c(-.02, .02)
  
  } 

if (meas == 'correlation'){
  if (suffix1 == '-same') mymeasures = c("correlation_trained_same", "correlation_untrained_same")
  if (suffix1 == '-different') mymeasures = c("correlation_trained_different", "correlation_untrained_different")
  if (suffix1 == '-all') mymeasures = c("correlation_trained_same", "correlation_untrained_same", "correlation_trained_different", "correlation_untrained_different", "correlation_trained_untrained")
  if (suffix1 == '-same') data$value = data$correlation_trained_same - data$correlation_untrained_same 
  else data$value = data$correlation_trained_different - data$correlation_untrained_different
  
#  data$value = data$correlation_trained_different - data$correlation_untrained_different 
  ylabel = "Correlation (z-scored)"
  ylimit = c(0, 0.2)
  ylimit_diff = c(-.005, .015)
} 

if (meas == 'mean_signal'){
  mymeasures = c("mean_signal_trained", "mean_signal_untrained")
  data$value = data$mean_signal_trained - data$mean_signal_untrained 
  ylabel = "% signal change"
  ylimit = c(-1, 3)
  ylimit_diff = c(-3, 3)
  
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

nrois = 4
# check data
plot(sort(table(data$subject))/nrois/2, las = 2) # how many sessions per subject
plot(sort(table(data$session))/nrois/2, las = 2) # how many sessions per timepoint

# relabel rois and separate hemispheres
control_labels = c("R_C1", "L_C1")
data = data %>% filter( !label %in% control_labels)
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
print(mymeasures)

data.melt = reshape2::melt(data, 
                 id.vars = c("subject", "session", "GROUP", "hemi", "label", 
                             "FD", "SYSTEM", "CONFIGURATION", "TRAINING"), #, "TRAINING.Q, TRAINING.A"
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
  theme(legend.title = element_blank(), legend.position = 'bottom',  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(myplot)


# hemi
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
  theme(text = element_text(size = 9),
        legend.title = element_blank(), legend.position = 'bottom', 
        legend.box = 'horizontal',
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

}

if (F) {
  WIDTH = 30; HEIGHT = 24; DPI = 1000
  meas = 'correlation'
  #meas = 'cosine'
  meas = 'euclidean'
  #meas = 'clf'
  suffix0 = 'mask-cross-perm'
  #suffix0 = 'mask-cross-derivatives'
  suffix0 = 'mask-cross'
  suffix1 = '-all'
  do_variability_analysis(meas, suffix0, suffix1)
}
