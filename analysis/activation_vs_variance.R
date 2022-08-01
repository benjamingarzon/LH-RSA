rm(list = ls()) 
library(reshape2)
library(ggplot2)
setwd("~/Software/LeftHand/analysis")
source("./plot_funcs.R")

incomplete_subjects = c("sub-lue5207", "sub-lue3203") # too few timepoints
suffix0 = 'mask-cross-runprew'
data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_', suffix0, '.csv')

data = read.table(data_file, header = T, sep = ',') %>% arrange(subject, session, hemi, label)
data$GROUP = "Intervention"
data$GROUP[grep("sub-lue.2", data$subject)] = "Control"

data = data %>% filter(!subject %in% incomplete_subjects)
control_labels = c("R_C1", "L_C1")
data = data %>% filter( !label %in% control_labels)
data$label  = gsub("R_", "Right ", data$label)
data$label  = gsub("L_", "Left ", data$label)
data$diff.xprod_trained = NA
data$diff.xprod_untrained = NA
data$diff.mean_signal_trained = NA
data$diff.mean_signal_untrained = NA
labels = unique(data$label)
subjects = sort(unique(data$subject))
mylist = NULL

get_slope = function(x, y){
  m = lm(y ~ x)
  return(coef(m)[2])
}


# correlation across sessions
for (mylabel in labels){
  for (mysubject in subjects){
    mydata = data %>% filter(subject == mysubject & label == mylabel)
    cor.untrained = cor.test(mydata$xproduct_grouped_untrained_same, mydata$mean_signal_untrained)$estimate
    cor.trained = cor.test(mydata$xproduct_grouped_trained_same, mydata$mean_signal_trained)$estimate
    sessions = sort(mydata$session)
    # print(sessions)
    if (1 %in% sessions) {

      diff.xprod_trained = mydata$xproduct_grouped_trained_same - mydata$xproduct_grouped_trained_same[1]
      data[data$subject == mysubject & data$label == mylabel,]$diff.xprod_trained = diff.xprod_trained

      diff.mean_signal_trained = mydata$mean_signal_trained - mydata$mean_signal_trained[1]
      data[data$subject == mysubject & data$label == mylabel,]$diff.mean_signal_trained = diff.mean_signal_trained
      
      diff.xprod_untrained = mydata$xproduct_grouped_untrained_same - mydata$xproduct_grouped_untrained_same[1]
      data[data$subject == mysubject & data$label == mylabel,]$diff.xprod_untrained = diff.xprod_untrained
      
      diff.mean_signal_untrained = mydata$mean_signal_untrained - mydata$mean_signal_untrained[1]
      data[data$subject == mysubject & data$label == mylabel,]$diff.mean_signal_untrained = diff.mean_signal_untrained

      # calculate slopes
      slope.xprod_trained = get_slope(mydata$session, mydata$xproduct_grouped_trained_same)
      slope.xprod_untrained = get_slope(mydata$session, mydata$xproduct_grouped_untrained_same)
      slope.mean_signal_trained = get_slope(mydata$session, mydata$mean_signal_trained)
      slope.mean_signal_untrained = get_slope(mydata$session, mydata$mean_signal_untrained)
      
    } else {
      slope.xprod_trained = slope.xprod_untrained = NA
      slope.mean_signal_trained = slope.mean_signal_untrained = NA
    }
    
    
    mylist = rbind(mylist, c(mysubject, unique(mydata$GROUP), mylabel, cor.trained, cor.untrained, 
                             slope.xprod_trained,
                             slope.xprod_untrained,
                             slope.mean_signal_trained,
                             slope.mean_signal_untrained))
  }
}

# correlation across labels of difference wrt to baseline
cor.diff.trained = NULL
cor.diff.untrained = NULL
for (mysession in sort(unique(data$session))[-1]){
#    for (mylabel in labels){
      cor.diff.trained.col = cor.diff.untrained.col = NULL
      for (mysubject in subjects){

        mydata = data %>% filter(session == mysession & subject == mysubject)
        indices = !is.na(mydata$diff.xprod_trained*mydata$diff.mean_signal_trained)
        cor.diff.trained.col = rbind(cor.diff.trained.col, 
                                     ifelse(sum(indices) > 0, cor.test(mydata$diff.xprod_trained[indices], 
                                              mydata$diff.mean_signal_trained[indices])$estimate, NA))
        indices = !is.na(mydata$diff.xprod_untrained*mydata$diff.mean_signal_untrained)
        
        cor.diff.untrained.col = rbind(cor.diff.untrained.col, 
                                      ifelse(sum(indices) > 0, cor.test(mydata$diff.xprod_untrained[indices], 
                                           mydata$diff.mean_signal_untrained[indices])$estimate, NA))
#      }
      }
  cor.diff.trained = cbind(cor.diff.trained, cor.diff.trained.col)
  cor.diff.untrained = cbind(cor.diff.untrained, cor.diff.untrained.col)
}

cor.diff.trained.melt = melt(cor.diff.trained) %>%filter(!is.na(value))
cor.diff.untrained.melt = melt(cor.diff.untrained) %>%filter(!is.na(value))
colnames(cor.diff.untrained.melt) = colnames(cor.diff.trained.melt) = c("subject", "session", "value")
cor.diff.trained.melt$measure = 'trained'
cor.diff.untrained.melt$measure = 'untrained'
cor.diff = rbind(cor.diff.trained.melt, cor.diff.untrained.melt)

cor.diff$session = cor.diff$session 
cor.diff$subject = subjects[cor.diff$subject]
cor.diff$group = "Intervention"
cor.diff$group[grep("sub-lue.2", cor.diff$subject)] = "Control"

cors.mean = cor.diff %>% group_by(group, measure, session)%>%
  summarise(val.mean = mean(value, na.rm = T), val.sem = sem(value)) 

# correlation across labels, separately for each group, averaging across subjects

myplot = ggplot(cors.mean, aes(x = session, y = val.mean,  ymin = val.mean - val.sem, group = measure, 
                               ymax = val.mean + val.sem, col = measure)) + 
  geom_point() + 
  geom_errorbar() + 
  geom_line() + 
  geom_hline(yintercept = 0, size = 0.3, linetype = 2) + 
  xlab('') + 
  ylab('Correlation') + 
  facet_grid('. ~ group') + 
  scale_colour_manual(values = myPalette2) +
  theme_lh + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=14, face="bold"),
        legend.box = 'vertical',
        legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "pt"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20, face="bold"),
        strip.text.x = element_text(size = 20, face="bold"),
        strip.background = element_blank())


print(myplot)

# correlation across sessions, separately for each group, averaging across subjects
cors = as.data.frame(mylist)[, 1:5]
colnames(cors) = c('subject', 'group', 'label', 'trained', 'untrained')
cors.melt = melt(cors, id.vars = c('subject', 'group', 'label'), variable.name = 'measure')
cors.melt$value = as.numeric(cors.melt$value)

cors.mean = cors.melt %>% group_by(group, label, measure)%>%
  summarise(val.mean = mean(value, na.rm = T), val.sem = sem(value)) 

myplot = ggplot(cors.mean, aes(x = measure, y = val.mean,  ymin = val.mean - val.sem, group = group, 
                               ymax = val.mean + val.sem, col = group)) + 
  geom_point() + 
  geom_errorbar() + 
  geom_line() + 
  geom_hline(yintercept = 0, size = 0.3, linetype = 2) + 
  xlab('') + 
  ylab('Correlation') + 
  scale_colour_manual(values = myPalette2) +
  facet_wrap('. ~ label', nrow = 2) +
  theme_lh + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=14, face="bold"),
        legend.box = 'vertical',
        legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "pt"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20, face="bold"),
        strip.text.x = element_text(size = 20, face="bold"),
        strip.background = element_blank())


print(myplot)

cors.test = cors.melt %>% group_by(group, label, measure) %>% summarise(pvalue = round(t.test(value)$p.value, digits = 3))
View(cors.test)

cors = as.data.frame(mylist)[, c(1:3, 6:9)]
colnames(cors) = c('subject', 'group', 'label', 'slope.xprod_trained', 'slope.xprod_untrained', 
                   'slope.mean_signal_trained', 'slope.mean_signal_untrained')

get_correlation = function(x, y) {
  x = as.numeric(x)
  y = as.numeric(y)
  valid = !is.na(x*y)
  x = x[valid]
  y = y[valid]
  rho = cor.test(x, y)$estimate
  return(rho)
}


cors.mean = cors %>% group_by(group, label)%>%
  summarise(trained = get_correlation(slope.xprod_trained, slope.mean_signal_trained),
            untrained = get_correlation(slope.xprod_untrained, slope.mean_signal_untrained))


cors.melt = melt(cors.mean, id.vars = c('group', 'label'), variable.name = 'measure')

myplot = ggplot(cors.melt, aes(x = measure, y = value,  group = group, col = group)) + 
  geom_point() + 
#  geom_errorbar() + 
  geom_line() + 
  geom_hline(yintercept = 0, size = 0.3, linetype = 2) + 
  xlab('') + 
  ylab('Correlation') + 
  scale_colour_manual(values = myPalette2) +
  facet_wrap('. ~ label', nrow = 2) +
  theme_lh + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=14, face="bold"),
        legend.box = 'vertical',
        legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "pt"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20, face="bold"),
        strip.text.x = element_text(size = 20, face="bold"),
        strip.background = element_blank())


print(myplot)

ff = subset(cors, group == "Intervention" & label == "Left PS")
plot(ff$slope.xprod_trained, ff$slope.mean_signal_trained)
plot(ff$slope.xprod_untrained, ff$slope.mean_signal_untrained)

