rm(list = ls())

sem = function(x) sd(x, na.rm = T)/sqrt(length(x))

data_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores.csv'
#output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/svm_acc.csv'
output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/alpha.csv'
data = read.table(data_file, header = T, sep = ',') %>% arrange(subject, session, hemi, label)

data$group = "Intervention"
data$group[grep("sub-lue.2", data$subject)] = "Control"

#data$value = data$alpha_trained - data$alpha_untrained #data$svm_acc
data$value = data$alpha_trained_PCM - data$alpha_untrained_PCM

data = data %>% filter(!is.infinite(value))

data.mean = data%>%
  group_by(session, group, hemi, label)%>%summarise(val.mean = mean(value, na.rm = T), 
                                       val.sem = sem(value))

myplot = ggplot(data.mean%>% filter(hemi == 'rh'), aes(
  x = session,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = group,
  group = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_wrap(~ label, ncol = 7) + 
  theme_classic()

print(myplot)

myplot = ggplot(data%>% group_by(session, group, hemi) %>% summarise(val.mean = mean(value),
                                                                    val.sem = sem(value)), aes(
  x = session,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = group,
  group = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  facet_grid(. ~ hemi) + 
  #ylim(0, 1) + 
  theme_classic()

print(myplot)

output.file = data%>%
  group_by(label)%>%
  summarise(value=mean(value)) %>% arrange(-value)

write.table(output.file, file = output_file, sep = ',', col.names = T, row.names = F)

# check data
plot(table(data$subject))
plot(table(data$session))

#alpha
plot(data$alpha_trained, data$alpha_untrained, pch = 20, cex = 0.2)
abline(0, 1)
alphas = dplyr::select(data, c(alpha_trained, alpha_untrained, alpha_trained_PCM, alpha_untrained_PCM))
alphas = alphas[rowMeans(alphas) > -Inf, ]

cor(alphas[complete.cases(alphas), ])
cor.test(data$alpha_trained, data$alpha_untrained)
cor.test(data$alpha_trained[!is.infinite(data$alpha_trained_PCM)], data$alpha_trained_PCM[!is.infinite(data$alpha_trained_PCM)])
cor.test(data$alpha_untrained[!is.infinite(data$alpha_untrained_PCM)], data$alpha_untrained_PCM[!is.infinite(data$alpha_untrained_PCM)])
cor.test(data$alpha_trained[!is.infinite(data$alpha_trained_PCM)], data$alpha_trained_PCM[!is.infinite(data$alpha_trained_PCM)])

plot(data$alpha_untrained, data$alpha_untrained_PCM)
