library(ggplot2)
library(dplyr)

rm(list = ls())
WD = '/data/lv0/MotorSkill/fmriprep/fmriprep/'
source('~/Software/ImageLMMR/ImageLMMR.R')
ANALYSIS_DIR = '/data/lv0/MotorSkill/fmriprep/analysis/'
TR = 1.2
paths = Sys.glob(file.path(
  WD,
  'sub-lue*',
  'ses-*',
  'func',
  '*_desc-confounds_regressors.tsv'
))

subject = session = run = FD = DVARS = adjR2 = NULL
for (i in seq(length(paths))) {
  path = paths[i]
  path_str = strsplit(path, '/')
  subject[i] = path_str[[1]][8]
  session[i] = as.numeric(strsplit(path_str[[1]][9], '-')[[1]][2])
  run[i] = as.numeric(substr(strsplit(path_str[[1]][11], '-')[[1]][5], 1, 1))
  pars = read.table(path, header = T)
  myFD = pars$framewise_displacement
  FD[i] = mean(as.numeric(myFD), na.rm = T)
  DVARS[i] = mean(as.numeric(pars$dvars), na.rm = T)
  print(paste(subject[i], session[i], run[i]))
  
  # see if the regressors predict FD
  X = NULL
  for (j in seq(6)){  
    TRs = (seq(length(myFD))-1)*TR
    myEV = TRs*0
    myfile = file.path(ANALYSIS_DIR, subject[i], sprintf('ses-%d/run%d/EV%d.csv', session[i], run[i], j))
    tryCatch({
        myonsets = read.table(myfile)
        myEV[sapply(myonsets$V1, function(x) which(TRs >= x)[1])] = 1
        X = cbind(X, myEV)
    }, error = function(e){})
  }
  if(!is.null(X)){
    colnames(X) = paste0("EV", seq(ncol(X)))
    X = as.data.frame(X)
    myFD[1] = 0
    X$y = as.numeric(myFD)
    model = lm(y ~ . , data = X)
    adjR2[i] = summary(model)$adj.r.squared
  } else adjR2[i] = NA
}

motion = data.frame(subject, session, run, FD, DVARS, adjR2)
motion$group = "Intervention"
motion$group[grep("sub-lue.2", motion$subject)] = "Control"

motion.orig = motion
write.table(motion.orig, file = '/data/lv0/MotorSkill/fmriprep/motion.csv')
# filter outliers
motion = motion %>% group_by(session) %>% mutate(FD = markoutliersIQR(FD)) %>% filter(!is.na(FD)) %>% ungroup()

sem = function(x) sd(x, na.rm = T)/sqrt(length(x))

motion.subject = motion%>%
  group_by(subject, session, group)%>%summarise(FD = mean(FD)) 

motion.mean = motion%>%
  group_by(subject, session, group)%>%summarise(FD = mean(FD), adjR2 = mean(adjR2)) %>% 
  group_by(session, group)%>%summarise(FD.mean = mean(FD, na.rm = T), FD.sem = sem(FD), FD.median = median(FD, na.rm = T),
                                       adjR2.mean = mean(adjR2, na.rm = T), adjR2.sem = sem(adjR2), adjR2.median = median(adjR2, na.rm = T))

myplot = ggplot(motion.subject, aes(
  x = session,
  y = FD,
  col = group,
  group = subject)) + 
  geom_line() + 
  geom_point() + 
  ylab('Mean Framewise Displacement') +
  xlab('Scanning session') +
  theme_classic()

print(myplot)

myplot = ggplot(motion.mean, aes(
  x = session,
  y = FD.mean,
  ymin = FD.mean - FD.sem, 
  ymax = FD.mean + FD.sem, 
  col = group,
  group = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  ylab('Mean Framewise Displacement') +
  xlab('Scanning session') +
#  ylim(c(0.1,0.25)) +
  theme_classic()

print(myplot)

myplot = ggplot(motion.mean, aes(
  x = session,
  y = adjR2.median,
  ymin = adjR2.median - adjR2.sem, 
  ymax = adjR2.median + adjR2.sem, 
  col = group,
  group = group
)) + geom_line() + 
  geom_point() + 
  geom_errorbar() + 
  ylab('Adj R2') +
  xlab('Scanning session') +
  #  ylim(c(0.1,0.25)) +
  theme_classic()

print(myplot)
