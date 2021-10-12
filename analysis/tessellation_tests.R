
rm(list = ls())

setwd("~/Software/LeftHand/analysis")
source('./plot_funcs.R')
source('./load_covariates.R')

suffix0 = 'tess-cross-runprew'
suffix1 = 'xnobis'
suffix = paste0(suffix0, suffix1)

data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_', suffix0, '.csv')
output_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/', suffix, '.csv')
figs_dir = '/data/lv0/MotorSkill/figs/tessellation'
data = read.table(data_file, header = T, sep = ',') %>% arrange(subject, session, hemi, label)

data$GROUP = "Intervention"
data$GROUP[grep("sub-lue.2", data$subject)] = "Control"

data$value = data$xnobis_trained_untrained

data = merge(data, covars.table %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
             by = c("subject", "session", "GROUP"))

data = left_join(data, motion %>% group_by(SUBJECT, TP) %>% 
                   summarise(FD = mean (FD)) %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
                 by = c("subject", "session"))

infos = NULL
par = "GROUPIntervention:TRAINING"
#par = "TRAINING"
mylabels = unique(data$label)
for (l in mylabels) {
  print(l)  
  model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*TRAINING + (1 |subject), data = data %>% 
                 filter(label == l))
  infos = rbind(infos, c(summary(model)$coefficients[par, ]))
}
par(mfrow = c(1, 2))
hist(infos[,'t value'], 30, xlab = 't')
hist(infos[,'Pr(>|t|)'], 30, xlab = 'p-value')
infos = as.data.frame(infos)
infos$label = mylabels
infos$value = infos[,'t value']
output.data = infos%>% arrange(-value) %>% select (label, value)

write.table(output.data, file = output_file, sep = ',', col.names = T, row.names = F)
