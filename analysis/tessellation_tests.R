
rm(list = ls())

library(mgcv)
k.SESS = 7
k.CONFIG = 7
smooth = "tp"
use_GAM = T

setwd("~/Software/LeftHand/analysis")
source('./plot_funcs.R')
source('./load_covariates.R')

suffix0 = 'tess-cross-runprew'
suffix1 = 'xnobis'
suffix = paste0(suffix0, suffix1)

data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_', suffix0, '.csv')
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
  if (!use_GAM) {
  model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*TRAINING + (1 |subject), data = data %>% 
                 filter(label == l))
  infos = rbind(infos, c(summary(model)$coefficients[par, ]))
  } else {
  model = gam(value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*TRAINING +
                s(subject, bs = "re") +
                s(TRAINING, by = interaction(GROUP), k = k.SESS, bs = smooth), 
              data = data %>% filter(label == l) %>% mutate(subject = factor(subject)), 
              method = 'REML')
  
  mysum = summary(model, freq = T)
  p.table = mysum$p.table
  p.table = cbind(p.table, df = df.residual(model))
  p.table = p.table[, c('Estimate', 'Std. Error', 'df', 't value', 'Pr(>|t|)')]
  infos = rbind(infos, p.table[par, ])
  }
  
}
par(mfrow = c(1, 2))
hist(infos[,'t value'], 30, xlab = 't')
hist(infos[,'Pr(>|t|)'], 30, xlab = 'p-value')
infos = as.data.frame(infos)
infos$label = mylabels
#infos$value = infos[,'t value']
#infos$value = p.adjust(infos[,'Pr(>|t|)'], method = 'fdr')
infos$value = 1 - p.adjust(infos[,'Pr(>|t|)'], method = 'fdr')
output.data = infos%>% dplyr::select (label, value)%>% arrange(-value) 

suffix = paste(suffix, par, sep = '-')
suffix = gsub(':', '_', suffix)
if (use_GAM) suffix = paste(suffix, 'GAM', sep = '-')
output_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/', suffix, '.csv')
write.table(output.data, file = output_file, sep = ',', col.names = T, row.names = F)
