rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)

# remove motion outliers...?

# difference alpha  trained > untrained increases with time and then goes to 0. More for intervention group
# check that not explained by noise!!
# use aroma to denoise??
# why is clf accuracy so high everywhere??
# correct for configuration x time
# correct for FD. FD vs type of 
# correct for number of correct trials
# do it for tessellation

setwd("~/Software/LeftHand/analysis")
source('./plot_funcs.R')
source('./load_covariates.R')

sem = function(x) sd(x, na.rm = T)/sqrt(length(x))
#data_file = '/data/lv0/MotorSkill/fmriprep/analysis/sub-lue5103/ses-4/surf/roi_scores_mask.csv'
#/surf/roi_scores_tess.csv' 
data_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_mask.csv'
#output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/clf_acc.csv'
output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/alpha_diff.csv'
#output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/crossnobis.csv'
data = read.table(data_file, header = T, sep = ',') %>% arrange(subject, session, hemi, label)

data$group = "Intervention"
data$group[grep("sub-lue.2", data$subject)] = "Control"
mymeasures = c("xnobis_trained_same", "xnobis_untrained_same", "xnobis_trained_different", "xnobis_untrained_different", "xnobis_trained_untrained")
#mymeasures = c("clf_acc_trained", "clf_acc_untrained")
#c("mean_signal_trained", "mean_signal_untrained") #c("alpha_trained", "alpha_untrained")
#mymeasures = c("crossnobis_trained", "crossnobis_untrained") #c("mean_signal_trained", "mean_signal_untrained") #c("alpha_trained", "alpha_untrained")

data$value = data$xnobis_untrained_same - data$xnobis_trained_same 
#data$value = data$mean_signal_trained
#data$value = data$clf_acc_trained - data$clf_acc_untrained
# data$value = data$clf_acc .25)/.25*100
#data$value = (data$clf_acc_perm - .25)/.25*100
#data$value = 0.5*(data$crossnobis_trained + data$crossnobis_untrained)
#data$value = data$crossnobis_trained - data$crossnobis_untrained
#data$value = data$alpha_trained_PCM - data$alpha_untrained_PCM
#data = subset(data, session != 1)
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

control_labels = c("R_C1", "L_C1")
incomplete_subjects = c("sub-lue5207", "sub-lue3202")
  
data = data %>% filter(!is.infinite(value))  %>% filter(!subject %in% incomplete_subjects, valid_runs >3) 

data = data %>% group_by(subject) %>% 
  mutate(value=markoutliersIQR(value)) %>% 
  filter(!is.na(value)) %>% ungroup()

data = merge(data, covars.table %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
             by = c("subject", "session"))

data = left_join(data, motion %>% group_by(SUBJECT, TP) %>% 
                   summarise(FD = mean (FD)) %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
                 by = c("subject", "session"))

data.melt = melt(data, #%>% filter( !label %in% control_labels), 
                 id.vars = c("subject", "session", "hemi", "label", "group", 
                             "GROUP", "FD", "SYSTEM", "CONFIGURATION", "TRAINING", "TRAINING.Q"),
                 variable.name = "MEASURE", 
                 value.name = "value") %>% mutate(SUBJECT = subject, value = as.numeric(value)) %>%
  filter(MEASURE %in% mymeasures) 


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
  group_by(session, group, hemi, label)%>%summarise(val.mean = mean(value, na.rm = T), 
                                       val.sem = sem(value))


myplot = ggplot(data.mean , aes( #%>% filter(hemi == 'rh')
  x = session,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = group,
  group = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
#  facet_wrap(CONFIGURATION ~ label, ncol = 5) + 
  facet_grid(hemi ~ label) + 
  theme_lh() + scale_colour_manual(values = myPalette)

print(myplot)

myplot = ggplot(data%>% filter( !label %in% control_labels) %>% group_by(session, group, hemi, CONFIGURATION) %>% 
                  summarise(val.mean = mean(value), val.sem = sem(value)), aes(
  x = session,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = group,
  group = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_grid(CONFIGURATION ~ hemi) + 
  #ylim(0, 1) + 
  theme_lh() + scale_colour_manual(values = myPalette)

print(myplot)

# hemi
data.mean = data.melt%>%
  group_by(session, group, label, MEASURE)%>%summarise(val.mean = mean(value.corr, na.rm = T), 
                                                    val.sem = sem(value.corr))%>% 
  mutate(group_measure = paste(group, MEASURE)) 

myplot = ggplot(data.mean, aes(
  x = session,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = MEASURE,
  group = group_measure,
  linetype = group
  )) + 
  geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_wrap(. ~label, ncol = 5) + 
  #facet_grid(. ~ hemi) + 
  theme_lh() + scale_colour_manual(values = myPalette) + theme(text = element_text(size = 8)) 

print(myplot)


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
