library(dplyr)
library(ggplot2)
library(reshape2)
library(viridis)
library(Rtsne)
library(plotly)
library(HardyWeinberg)
library(gg3D)

do_session_variability_analysis = function(mymetric, suffix0, suffix1, mylabels){
suffix = paste0(suffix0, suffix1)

figs_dir = '/data/lv0/MotorSkill/figs/session_variability'

if (suffix1 == '-same') seqs_train = c('trained_same', 'untrained_same')
if (suffix1 == '-different') seqs_train = c('trained_different', 'untrained_different', 'trained_untrained')
if (suffix1 == '-all') seqs_train = c('trained_same', 'untrained_same', 'trained_different', 'untrained_different', 'trained_untrained')

ylabels = c('Correlation (z-scored)', 'Cosine (z-scored)', 'Cross-nobis distance', 'Squared euclidean distance')
names(ylabels) = c('xcorrelation', 'xcosine', 'xnobis', 'xeuclidean')

########################
# Prepare data
########################
data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/across_session_scores_', suffix0, '.csv')
data = read.table(data_file, header = T, sep = ',')

data$group = "Intervention"
data$group[grep("sub-lue.2", data$subject)] = "Control"
data = data[complete.cases(data), ]

data = data %>% group_by(subject) %>% 
  mutate(value=markoutliersIQR(value)) %>% 
  filter(!is.na(value)) %>% ungroup()

data = data %>% mutate(session_diff = abs(session_train - session_test))

# remove first wave
data = data[ -grep('sub-lue1', data$subject), ]

# relabel
data$label  = gsub("R_", "Right ", data$label)
data$label  = gsub("L_", "Left ", data$label)
data$label  = gsub("C1", "Control Region", data$label)
#data$uni_label  = gsub("Right ", "Average ", data$label)
#data$uni_label  = gsub("Left ", "Average ", data$uni_label)

data.mean = data%>%
  group_by(group, label, session_train, session_test, session_diff, metric, seq_train)%>% 
  mutate(value = 1 - ifisherz(value)) %>%
  dplyr::summarise(val.mean = mean(value, na.rm = T), val.sem = sem(value))%>% 
  mutate(train_group = paste(seq_train, group))

data.subset.session = data.mean %>% filter(#session_test != session_train &
                                     metric == mymetric & 
                                     label %in% mylabels & 
                                     seq_train %in% seqs_train) %>%
  mutate(seq_train = clean_measure_names(seq_train, keep_first = T))

########################
# Plots
########################
myplot.session = ggplot(data.subset.session, aes( 
  x = session_train,
  y = val.mean,
  ymin = val.mean - val.sem, 
  ymax = val.mean + val.sem, 
  col = seq_train,
  group = train_group,
  linetype = group
)) + geom_line(size = 0.5) + 
  geom_point(size = 0.5) + 
  geom_errorbar() + 
  facet_grid( label ~ session_test) + 
  theme_lh() + theme(legend.position = 'bottom') +
  ylab(ylabels[mymetric]) +
  xlab('Test session') + 
  scale_colour_manual(values = myPalette) +
  theme(legend.title = element_blank(), legend.position = 'bottom',  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(myplot.session)


myplot = ggplot(data.subset.session %>% filter(session_diff <= 3), aes( 
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
  facet_grid( label + group ~ session_train) + 
  theme_lh() + theme(legend.position = 'bottom') + scale_colour_manual(values = myPalette)

print(myplot)


#data.subset = data.mean %>% filter(metric == mymetric & 
#                                     seq_train %in% seqs_train)

myplot.matrix = ggplot(data.subset.session, aes( 
  x = session_test,
  y = session_train,
  fill = val.mean)) +  
  geom_tile() + 
  facet_grid(label ~ seq_train + group) + 
  theme_lh() + theme(legend.position = 'bottom') + scale_fill_viridis() + 
  xlab('Test session') + 
  ylab('Test session') +
  labs(fill = ylabels[mymetric])

print(myplot.matrix)

ggsave(
  file.path(figs_dir, paste(mymetric, suffix, 'session.png', sep = '-')),
  plot = myplot.session,
  width = WIDTH,
  height = HEIGHT,
  dpi = DPI,
  units = 'cm'
)

ggsave(
  file.path(figs_dir, paste(mymetric, suffix, 'matrix.png', sep = '-')),
  plot = myplot.matrix,
  width = WIDTH,
  height = HEIGHT,
  dpi = DPI,
  units = 'cm'
)

# Some statistics
# in intervention, trained decorrelate fmore the farther they are
data.subset = data %>% filter(session_test != session_train & 
#                                session_train == 1 &
                                metric == mymetric & 
                                label == 'Left PS' & 
                                seq_train %in% c("trained_different", "untrained_different") )

# are the matrices symetric?
if (suffix1 == '-all') {
model = lmer(value ~ scale(session_diff)*seq_train*group + (1| subject), data = data.subset)
summary(model)

data.cmd.all = NULL
#for (mylabel in unique(data.mean$label)){
#for (mylabel in c('Right SPL', 'Right PM', 'Right PS')){
for (mylabel in mylabels){
    
# mds - only with suffix1  == all
data.mds = data.mean %>% filter(metric == mymetric & 
                                  seq_train %in% seqs_train & 
                                  label == mylabel) %>% 
  arrange(label, group, seq_train, session_train, session_test)
data.mds.2 = data.mds %>% arrange(label, group, seq_train, session_test, session_train)
data.mds$val.mean.train = data.mds$val.mean
data.mds$val.mean.test = data.mds.2$val.mean
data.mds$val.mean = 0.5*(data.mds$val.mean.train + data.mds$val.mean.test)

mysessions = seq(7)
NSESSIONS = length(mysessions)
sim.mat = matrix(0, NSESSIONS*6, NSESSIONS*6)
for (mygroup in c('Control', 'Intervention')) {
  for (sess_train in seq(NSESSIONS))
  {
    for (sess_test in seq(NSESSIONS)){
      mydata = subset(data.mds, session_train == mysessions[sess_train] & session_test == mysessions[sess_test])
      block = matrix(0, 6, 6)
      
      trained_same = mydata$val.mean[mydata$seq_train == 'trained_same' & mydata$group == mygroup ]
      trained_different = mydata$val.mean[mydata$seq_train == 'trained_different' & mydata$group == mygroup ]
      untrained_same = mydata$val.mean[mydata$seq_train == 'untrained_same' & mydata$group == mygroup ]
      untrained_different = mydata$val.mean[mydata$seq_train == 'untrained_different' & mydata$group == mygroup ]
      trained_untrained = mydata$val.mean[mydata$seq_train == 'trained_untrained' & mydata$group == mygroup ]

      if (sess_test != sess_train) {
        block[1, 1] = block[2, 2] = block[5, 5] = trained_same
        if (length(untrained_same)>0){
          block[3, 3] = block[4, 4] = block[6, 6] = untrained_same
        } else {
          block[3, 3] = block[4, 4] = block[6, 6] = untrained_different
        }
      } 
      block[1, 2] = trained_same

      block[1:2, 3] = trained_untrained

      block[1:2, 4] = trained_untrained
      
      if (length(untrained_same)>0){
        block[3, 4] = untrained_same
      } else {
        block[3, 4] = untrained_different
      }

      block[1:2, 5] = trained_different
      block[3:4, 5] = trained_untrained
      
      block[1:2, 6] = trained_untrained
      block[3:4, 6] = untrained_different
      block[5, 6] = trained_untrained
      
      block = block + t(block)
      seq_types = c('trained_A1', 'trained_A2', 'untrained_C1', 'untrained_C2', 'trained_B', 'untrained_D')
      sim.mat[((sess_train - 1)*6 + 1) : ((sess_train - 1)*6 + 6),
              ((sess_test - 1)*6 + 1) : ((sess_test - 1)*6 + 6)] = block
    }
  }
  dist.mat = 1 - ifisherz(sim.mat)
  dist.mat[eye(nrow(dist.mat)) == 1] = 0
  ncoord = 2
  cmd = cmdscale(dist.mat, k = ncoord, eig = T)
  data.cmd = as.data.frame(cmd$points)
  colnames(data.cmd) = c('x', 'y', 'z')[1:ncoord]
  data.cmd$seq_types = rep(seq_types, NSESSIONS)
  data.cmd$session = as.vector(sapply(mysessions, function(x) rep(x, 6)))

  data.cmd$label = mylabel
  data.cmd$group = mygroup
  
  data.cmd.all = rbind(data.cmd.all, data.cmd)
} # group
} #label


# reference everyone to origin
if (T)
{for (group in unique(data.cmd.all$group)){
 for (label in unique(data.cmd.all$label)){
    for (seq_type in unique(data.cmd.all$seq_types)) {
    indices = data.cmd.all$seq_types == seq_type & data.cmd.all$label == label & data.cmd.all$group == group
    data.cmd.all[ indices, c('x', 'y', 'z')[1:ncoord]] =
      data.cmd.all[ indices, c('x', 'y', 'z')[1:ncoord]] - 
      matrix(rep(unlist(data.cmd.all[indices & data.cmd.all$session == 1, c('x', 'y', 'z')[1:ncoord] ]), NSESSIONS), 
             ncol = ncoord, byrow = T)
  }
  }
}
}
data.cmd.all$seq_types[data.cmd.all$seq_types == 'untrained_D' & data.cmd.all$session %in% c(2, 5)] = 'untrained_E'
data.cmd.all$seq_types[data.cmd.all$seq_types == 'untrained_D' & data.cmd.all$session %in% c(3, 6)] = 'untrained_F'

View(data.cmd.all)

myplot.prog = ggplot(data.cmd.all, aes(x = x, y = y, col = session, group = seq_types, shape = seq_types)) + 
 # axes_3D() + 
  geom_path() + 
  geom_point(size = 4) + 
#  stat_3D(size = 4) + 
  xlab("Coordinate 1") + 
  ylab("Coordinate 2") + 
#  zlab("Coordinate 3") + 
  theme_lh() +
  scale_colour_viridis() + 
  scale_shape_manual(values=1:8) + 
  facet_grid(group ~ label)
#facet_wrap(group ~ label, scales = "free", nrow = 2)

ggsave(
  file.path(figs_dir, paste(mymetric, suffix, 'progression.png', sep = '-')),
  plot = myplot.prog,
  width = WIDTH,
  height = HEIGHT,
  dpi = DPI,
  units = 'cm'
)
}

}

if (T) {
  ########################
  # Options
  ########################
  WIDTH = 30
  HEIGHT = 24
  DPI = 1000
  mymetric = 'xcorrelation'
  suffix0 = 'mask-cross'
  suffix1 = '-all'
  mylabels = c('Right SPL', 'Right PS', 'Right PM', 'Left PM', 'Left SPL') #, 'Right Control Region')
  do_session_variability_analysis(mymetric, suffix0, suffix1, mylabels)
}
