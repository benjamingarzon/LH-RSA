rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(viridis)
library(Rtsne)

setwd("~/Software/LeftHand/analysis")
source('./plot_funcs.R')
source('./load_covariates.R')

sem = function(x) sd(x, na.rm = T)/sqrt(length(x))
data_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/across_session_scores_mask-cross.csv'
data = read.table(data_file, header = T, sep = ',')

data$group = "Intervention"
data$group[grep("sub-lue.2", data$subject)] = "Control"
data = data[complete.cases(data), ]

data = data %>% group_by(subject) %>% 
  mutate(value=markoutliersIQR(value)) %>% 
  filter(!is.na(value)) %>% ungroup()

data = data %>% mutate(session_diff = abs(session_train - session_test))
data.mean = data%>%
  group_by(group, label, session_train, session_test, session_diff, metric, seq_train)%>%dplyr::summarise(val.mean = mean(value, na.rm = T), 
                                                    val.sem = sem(value))%>% mutate(train_group = paste(seq_train, group))

mymetric = 'xcorrelation'
mylabels = c('R_SPL', 'R_C1') #'R_C1',
#seqs_train = c('same', 'different')
#seqs_train = c('trained_different', 'untrained_different') #, 'trained_untrained')
seqs_train = c('trained_same', 'untrained_same') #, 'trained_untrained')

data.subset = data.mean %>% filter(session_test != session_train &
                                     metric == mymetric & 
                                     label %in% mylabels & 
                                     seq_train %in% seqs_train) 

myplot = ggplot(data.subset, aes( 
  x = session_train,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = seq_train,
  group = train_group,
  linetype = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_grid( group + label ~ session_test) + 
  theme_lh() + theme(legend.position = 'bottom') + scale_colour_manual(values = myPalette)

print(myplot)


myplot = ggplot(data.subset %>% filter(session_diff <= 3), aes( 
  x = session_diff,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = seq_train,
  group = train_group,
  linetype = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_grid( group + label ~ session_train) + 
  theme_lh() + theme(legend.position = 'bottom') + scale_colour_manual(values = myPalette)

print(myplot)


#data.subset = data.mean %>% filter(metric == mymetric & 
#                                     seq_train %in% seqs_train)

myplot = ggplot(data.subset, aes( 
  x = session_test,
  y = session_train,
  fill = val.mean)) +  
  geom_tile() + 
  facet_grid(label + metric ~ seq_train + group) + 
  theme_lh() + theme(legend.position = 'bottom') + scale_fill_viridis()

print(myplot)

# some stats
data.subset = data %>% filter(session_test != session_train & 
#                                session_train == 1 &
                                metric == mymetric & 
                                label == 'R_C1' & 
                                seq_train %in% seqs_train) 

model = lmer(value ~ session_diff*seq_train*group + (1| subject), data = data.subset)
summary(model)
