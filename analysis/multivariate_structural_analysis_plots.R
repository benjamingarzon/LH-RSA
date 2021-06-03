rm(list = ls())

######################################
# Plot results of multivariate analyses
######################################
library(freesurfer)
library(dplyr)
library(ggplot2)
source("./structural_analysis_funcs.R")
source("./plot_funcs.R")

#load behavioural parameters
load('~/Data/LeftHand/Lund1/responses/MT_learning_parameters.RData')

# redo label
alphalevel = 0.05
DPI = 1000
FEAT.THR = 0.8 # prop of coefs of equal sign
ISTHICKNESS = T
MYYLABEL = ifelse(ISTHICKNESS, "% increase in cortical thickness", "% increase in grey matter volume")
USEVALIDFEATURES = F
OUTPUTDIR='multivariate-tess'
inflated.rh = freesurfer_read_surf('/usr/local/freesurfer/subjects/fsaverage/surf/rh.inflated')
inflated.lh = freesurfer_read_surf('/usr/local/freesurfer/subjects/fsaverage/surf/lh.inflated')


if (ISTHICKNESS) {
  FIGS_DIR='/home/benjamin.garzon/Data/LeftHand/Lund1/figs/structure/multivar-thickness'
  
  VOXELORVERTEX = 'Vertex'
  DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results_unsmoothed'
  load(file = file.path(DATADIR, 'tests', OUTPUTDIR, 'results.rds'))
  
  # collect and plot results
  accuracy.rh = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'rh.accuracy'), header = T, sep = ',', na.strings = "NA", strip.white = T)
  accuracy.lh = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'lh.accuracy'), header = T, sep = ',', na.strings = "NA", strip.white = T)
  accuracy.rh$hemi = "rh"
  accuracy.lh$hemi = "lh"
  
  accuracy = rbind(accuracy.rh, accuracy.lh)
  
  # compare permutation binomial distribution
  accuracy.rh.binom = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'rh.accuracy.binom'), header = T, sep = ',', na.strings = "NA", strip.white = T)
  accuracy.lh.binom = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'lh.accuracy.binom'), header = T, sep = ',', na.strings = "NA", strip.white = T)
  accuracy.rh.binom$hemi = 'rh'
  accuracy.lh.binom$hemi = 'lh'
  accuracy.binom = rbind(accuracy.rh.binom, accuracy.lh.binom)
  
  train_res_list = c(results.multivariate.rh$train_res_list, results.multivariate.lh$train_res_list)
  
} else {
  # load results of VBM
  FIGS_DIR='/home/benjamin.garzon/Data/LeftHand/Lund1/figs/structure/multivar-VBM'
  VOXELORVERTEX = 'Voxel'
  
  DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/'
  load(file = file.path(DATADIR, 'tests', OUTPUTDIR, 'results.rds'))
  train_res_list = c(results.multivariate.vol$train_res_list)
  
  accuracy = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'MNI.accuracy'), header = T, sep = ',', na.strings = "NA", strip.white = T)
  accuracy$hemi = sapply(as.character(accuracy$label), function(x) strsplit(x, '\\.')[[1]][1])
  accuracy$label = sapply(as.character(accuracy$label), function(x) strsplit(x, '\\.')[[1]][2])
  
  # compare permutation binomial distribution
  accuracy.binom = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'MNI.accuracy.binom'), header = T, sep = ',', na.strings = "NA", strip.white = T)
  accuracy.binom$hemi = accuracy$hemi 
  accuracy.binom$label = accuracy$label
  
}

######################################
# Compute p-values
######################################
means = colMeans(dplyr::select(accuracy, -c(label, hemi)))
maxs = apply(dplyr::select(accuracy, -c(label, hemi)), 2, max)

pvals = apply(dplyr::select(accuracy, -c(label, hemi)), 1, function(x) mean(x[1] <= x))
pvals.corr = sapply(accuracy$value, function (x) mean(x <= maxs)) 
res = cbind(dplyr::select(accuracy, c(label, hemi, value)), pvals = pvals, pvals.corr = pvals.corr )%>% arrange(pvals.corr, pvals, desc(accuracy))
#print(means)
#print(maxs)
res.significant = subset(res, pvals <= alphalevel)
res.corr = subset(res, pvals.corr <= alphalevel)


which.res = res.corr
data.max = expand.grid(label = accuracy$label, m = maxs[-1]) %>% filter(label %in% which.res$label)
data.max = merge(data.max, res.corr, by = 'label')
data.max$full_label = paste(data.max$hemi, data.max$label, '.')
which.res$full_label = paste(which.res$hemi, which.res$label, '.')

myplot = ggplot() + 
  geom_violin(data = data.max, aes (x = reorder(full_label, -value), y = m), alpha = 0.5) +  
  geom_point(data = which.res, aes(x = reorder(full_label, -value), y = value, col = hemi), size = 2) + 
  ylim(0.5, 0.7) + 
  scale_colour_manual(values = myPalette) + theme_lh + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') + 
  xlab('Region') + ylab('Accuracy') #+ facet_grid(. ~ hemi)
print(myplot)
ggsave(file.path(FIGS_DIR, 'Accuracy.png'), dpi = DPI)

write.table(filter(res.corr, hemi == "rh") %>% dplyr::select(c(label, value)), file.path(DATADIR, 'tests', OUTPUTDIR, 'rh.accuracy.results.csv'), 
            col.names = T, sep = ',', row.names = F)
write.table(filter(res.corr, hemi == "lh") %>% dplyr::select(c(label, value)), file.path(DATADIR, 'tests', OUTPUTDIR, 'lh.accuracy.results.csv'), 
            col.names = T, sep = ',', row.names = F)
write.table(dplyr::select(res.corr, c(hemi, label, value)), file.path(DATADIR, 'tests', OUTPUTDIR, 'accuracy.csv'), 
            col.names = T, sep = ',', row.names = F)

maxs.binom = apply(dplyr::select(accuracy.binom, -c(label, hemi)), 2, max)

breaks = seq(0.4, 0.7, 0.01)
hist(maxs[-1], breaks = breaks, col = 'red', ylim = c(0, 15)) # red
hist(maxs.binom[-1], breaks = breaks, add = T, col = rgb(0, 1, 0, 0.4)) # green


######################################
# Plot differences
######################################

## original data, subject averages
## pca data, subject averages and 
# how to compare different regions?
# compare stability of the features and weight
# try with Glasser parcellation...

#####
# change the coefs so that they can be interpreted
#####

# check distribution for different models
which.model = sapply(train_res_list, function(x) x$prepared_data$which.model)
myperms = sapply(train_res_list, function(x) x$perm) == 0
mylabels = as.character(sapply(train_res_list, function(x) x$label))[myperms]

if (ISTHICKNESS) mylabels = paste(paste(accuracy$hemi, mylabels, sep = '.'))

mypcadim = sapply(train_res_list, function(x) ncol(x$pca$x))[myperms]

which.model = matrix(unlist(which.model[myperms]), ncol = length(mylabels))
colnames(which.model) = mylabels
mygroup = train_res_list[[1]]$y
model_types = c('c0', 'l0', 'l1', 'q1', 'q2')
model.prop.control = apply(which.model[mygroup == 1, ], 2, function(x) table(factor(x, levels = model_types))) #/nrow(which.model)
model.prop.exp = apply(which.model[mygroup == 2, ], 2, function(x) table(factor(x, levels = model_types))) #/nrow(which.model)
model.prop = apply(which.model, 2, function(x) table(factor(x, levels = model_types)))/nrow(which.model)


par(mfrow = c(2,3))
for (model_type in model_types)
{ 
  plot(model.prop.control[model_type, ], type = 'l', main = model_type)
  lines(model.prop.exp[model_type, ], col = "red")
  #  plot(model.prop.control[model_type, ],  
  #       model.prop.exp[model_type, ], main = model_type, type = 'p')
  
}


# count most stable features
myvalidfeatures = sapply(train_res_list[myperms], function(z) apply(sign(z$coefs.iter), 2, function(x) {y = x[x!=0]; max(table(y)/length(x), na.rm = T) > FEAT.THR}))
nfeatures = sapply(myvalidfeatures, function(x) sum(x))
par(mfrow = c(1, 2))
hist(nfeatures, 20)
plot(nfeatures, mypcadim)
names(myvalidfeatures) = mylabels 


data.proj.all = NULL
for (j in seq(nrow(res.corr))){
  label = res.corr$label[j]
  hemi = res.corr$hemi[j]
  if (ISTHICKNESS) {
    if (hemi == "rh") {
      train_res_list = results.multivariate.rh$train_res_list
      coords = inflated.rh$vertices
    } else {
      train_res_list = results.multivariate.lh$train_res_list
      coords = inflated.lh$vertices
    }
  } else {
    train_res_list = results.multivariate.vol$train_res_list
    label = paste(hemi, label, sep = '.')
    if (hemi == "rh") {
      coords = inflated.rh$vertices
    } else {
      coords = inflated.lh$vertices
    }
  }
  
  print(label)
  fulllabel = ifelse(ISTHICKNESS, paste(hemi, label, sep = '.'), label)
  # take not permutated instance (perm == 0)
  w = which(sapply(train_res_list, function(x) x$label == label & x$perm == 0))
  train_res = train_res_list[[w]]
  mask.indices = which(train_res$mask)
  valid.features = unlist(myvalidfeatures[fulllabel])
  scores = train_res$pca$x[, valid.features, drop = F]
  loadings = train_res$pca$rotation[, valid.features, drop = F]
  
  c1 = scores[train_res$y == 1,]
  c2 = scores[train_res$y == 2,] # experimental
  
  if(USEVALIDFEATURES) {
    # show results only considering features that separate cases
    x1 = loadings %*% t(c1)
    x2 = loadings %*% t(c2)
  } else {
    x1 = t(train_res$prepared_data$mycoefs[train_res$y == 1,])
    x2 = t(train_res$prepared_data$mycoefs[train_res$y == 2,])
  }
  
  m1 = rowMeans(x1)
  m2 = rowMeans(x2)
  s1 = apply(x1, 1, se)
  s2 = apply(x2, 1, se)
  
  # unfold coefs
  indices = t(sapply(rownames(loadings), function(x)
    unlist(strsplit(x, '_'))))
  colnames(indices) = c("ORDER", "VERTEX")
  data.proj = as.data.frame(cbind(
    indices,
    m1 = m1,
    m2 = m2,
    s1 = s1,
    s2 = s2
  )) %>% mutate_all(as.character) %>% mutate_all(as.numeric)
  
  train_res$c1 = c1
  train_res$c2 = c2
  train_res$s1 = s1
  train_res$s2 = s2
  train_res$data.proj = data.proj
  
  data.proj$INDEX = mask.indices[data.proj$VERTEX] # index within mask
  data.proj$x = coords[data.proj$INDEX, 1] 
  data.proj$y = coords[data.proj$INDEX, 2]
  data.proj$z = coords[data.proj$INDEX, 3]
  data.proj$LABEL = label
  data.proj$HEMI = hemi  
  #visualize
  myids = c("VERTEX", "ORDER", "INDEX", "x", "y", "z", "LABEL", "HEMI")
  data.proj.melt.m = reshape2::melt(data.proj, id.vars = myids, 
                                    measure.vars = c("m1", "m2"), variable.name = "GROUP", value.name = "m") %>% 
    mutate(GROUP = ifelse(GROUP == "m1", "Control", "Intervention")) %>% arrange(LABEL, HEMI, GROUP, VERTEX, ORDER, INDEX)
  g1 = data.proj.melt.m$GROUP == "Intervention" & data.proj.melt.m$ORDER == 1
  g2 = data.proj.melt.m$GROUP == "Control" & data.proj.melt.m$ORDER == 1
  g3 = data.proj.melt.m$GROUP == "Intervention" & data.proj.melt.m$ORDER == 2
  g4 = data.proj.melt.m$GROUP == "Control" & data.proj.melt.m$ORDER == 2
  
  data.proj.melt.m$rank[g1] = 
    rank(data.proj.melt.m$m[g1])
  data.proj.melt.m$rank[g2] = data.proj.melt.m$rank[g1]
  data.proj.melt.m$rank[g3] = 
    rank(data.proj.melt.m$m[g3])
  data.proj.melt.m$rank[g4] = data.proj.melt.m$rank[g3]
  
  data.proj.melt.s = reshape2::melt(data.proj, id.vars = myids, 
                                    measure.vars = c("s1", "s2"), variable.name = "GROUP", value.name = "s")%>% 
    mutate(GROUP = ifelse(GROUP == "s1", "Control", "Intervention"))
  
  data.proj.melt = merge(data.proj.melt.m, data.proj.melt.s, by = c(myids, "GROUP")) 
  data.proj.all = rbind(data.proj.all, data.proj.melt)
  
}

if (ISTHICKNESS) {
  data.proj.all = merge(data.proj.all, dplyr::select(accuracy, c(hemi, label, value))%>% dplyr::rename(ACCURACY = value, HEMI = hemi, LABEL = label), by = c("HEMI", "LABEL")) 
} else {
  data.proj.all = merge(data.proj.all, dplyr::select(accuracy, c(hemi, label, value))%>% dplyr::rename(ACCURACY = value) %>% 
                          mutate(LABEL = paste(hemi, label, sep = '.')), by = c("LABEL")) 
}
data.proj.wide = reshape(data.proj.all, idvar = c("VERTEX", "LABEL", "HEMI", "GROUP"), timevar = "ORDER", direction = "wide")
data.proj.all$ORDER = ifelse(data.proj.all$ORDER == 1, "Linear term", "Quadratic term")

if (ISTHICKNESS) {
  data.proj.all$FULL_LABEL = paste(data.proj.all$HEMI, data.proj.all$LABEL, sep ='.')
  data.proj.wide$FULL_LABEL = paste(data.proj.wide$HEMI, data.proj.wide$LABEL, sep ='.')
} else {
  data.proj.all$FULL_LABEL = data.proj.all$LABEL
  data.proj.wide$FULL_LABEL = data.proj.wide$LABEL
}
######################################
# Plot differences
######################################

# visualize patches 
myplot = ggplot(data.proj.all, aes ( x = y, y = z, col = m)) + geom_point() + 
  scale_colour_viridis_c() + facet_grid(FULL_LABEL + GROUP ~ ORDER)
print(myplot)

# group vertices together
myplot = ggplot(data.proj.all, aes ( x = FULL_LABEL, 
                                     y = m, 
                                     ymin = m - s, 
                                     ymax = m + s,
                                     fill = as.factor(GROUP))) + 
  #geom_line(size = 0.1, alpha = 0.1) + 
  geom_violin() +
  xlab('Region') +
  ylab('Coefficient') + 
  facet_grid(as.factor(ORDER) ~ .) + 
  scale_fill_manual(values = myPalette) + theme_lh +  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.title = element_blank(), legend.position = "bottom") 

print(myplot)

ggsave(file.path(FIGS_DIR, 'CoefsvsRegions.png'), dpi = DPI)

# show by vertex/voxel
myplot = ggplot(data.proj.all, aes ( x = rank, 
                                     y = m, 
                                     ymin = m - s, 
                                     ymax = m + s, 
                                     col = GROUP, 
                                     fill = GROUP,
                                     group = GROUP)) + 
  geom_line(size = 0.5) + 
  geom_ribbon(alpha = 0.5, size = 0) +
  facet_grid( as.factor(ORDER) ~ reorder(FULL_LABEL, -ACCURACY)) + xlab(VOXELORVERTEX) + ylab('Coefficient') + 
  scale_colour_manual(values = myPalette) + theme_lh + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size=1), legend.title = element_blank(), legend.position = "bottom")
print(myplot)

ggsave(file.path(FIGS_DIR, 'IndividualRegions.png'), dpi = DPI)

# scatterplots
myplot = ggplot(data.proj.wide, aes ( x = m.1, 
                                      y = m.2,  
                                      col = GROUP)) + 
  geom_point(size = 0.5) + 
  facet_wrap(FULL_LABEL ~ .) + xlab('Linear term') + ylab('Quadratic term') + 
  scale_colour_manual(values = myPalette) + theme_lh + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.title = element_blank(), legend.position = "bottom") 
print(myplot)

ggsave(file.path(FIGS_DIR, 'ScatterplotLinearvsQuadratic.png'), dpi = DPI)

# now everyone in the same plot
myplot = ggplot(data.proj.wide, aes ( x = m.1, 
                                      y = m.2,  
                                      col = FULL_LABEL)) + 
  geom_point(size = 0.5) + 
  xlab('Linear term') + ylab('Quadratic term') + 
  facet_grid(GROUP ~ .) + 
  theme_lh + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.title = element_blank(), legend.position = "bottom") 
print(myplot)

ggsave(file.path(FIGS_DIR, 'ScatterplotLinearvsQuadraticTogether.png'), dpi = DPI)

# order 1 vs order 2
myplot = ggplot(data.proj.all, aes ( x = rank, 
                                     y = m, 
                                     ymin = m - s, 
                                     ymax = m + s, 
                                     col = GROUP, 
                                     fill = GROUP,
                                     group = GROUP)) + 
  geom_line(size = 0.5) + 
  geom_ribbon(alpha = 0.5, size = 0) +
  facet_grid(as.factor(ORDER)~ HEMI) + xlab(VOXELORVERTEX) + ylab('Coefficient') + 
  scale_colour_manual(values = myPalette) + theme_lh + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) 
print(myplot)


######################################
# Associations with behavioural learning curves and timecourse plots
######################################
#c("r0_norm", "e0_norm", "p0_norm", "kr_norm", "ke_norm", "kp_norm")
parnames = c("kr_norm", "ke_norm", "kp_norm")
par(mfrow = c(ceiling(nrow(res.corr)/2), length(parnames)))
XX.all = NULL
myplots = NULL
j = 1
for (j in seq(nrow(res.corr))){
  label = res.corr$label[j]
  hemi = res.corr$hemi[j]
  if (ISTHICKNESS) {
    if (hemi == "rh") {
      train_res_list = results.multivariate.rh$train_res_list
      coords = inflated.rh$vertices
    } else {
      train_res_list = results.multivariate.lh$train_res_list
      coords = inflated.lh$vertices
    }
  } else {
    train_res_list = results.multivariate.vol$train_res_list
    label = paste(hemi, label, sep = '.')
    if (hemi == "rh") {
      coords = inflated.rh$vertices
    } else {
      coords = inflated.lh$vertices
    }
  }

  # take not permutated instance (perm == 0)
  w = which(sapply(train_res_list, function(x) x$label == label & x$perm == 0))
  train_res = train_res_list[[w]]
  fulllabel = ifelse(ISTHICKNESS, paste(hemi, label, sep = '.'), label)
  valid.features = unlist(myvalidfeatures[fulllabel])
  scores = train_res$pca$x[, valid.features, drop = F]
  #loadings = train_res$pca$rotation[, valid.features, drop = F]
  
  ###################################### remove
  myvars = c("SUBJECT",
             "TP",
             "FD", 
             "SYSTEM",
             "TRAINING",
             "TRAINING.Q",
             "TRAINING.C",
             "GROUP.NUM", 
             "GROUP")
  XX = cbind(train_res$prepared_data$DATA[myvars], train_res$prepared_data$X) 
  XX.melt = reshape::melt(XX, id.vars = myvars, variable_name = "VERTEX")
  XX.init = unique(XX.melt[c("SUBJECT", "GROUP.NUM")])
  XX.melt$FULL_LABEL = fulllabel
  XX.mean = XX.melt %>% mutate(WAVE = substr(SUBJECT, 4, 4)) %>%group_by(SUBJECT, WAVE, TP, FD, SYSTEM, TRAINING, TRAINING.Q, GROUP, FULL_LABEL) %>% 
    summarise(y = mean(value)) # average across vertices
  myplots[[j]] = plot_data(XX.mean, MYYLABEL)
  j = j + 1
  # XX.all = rbind(XX.all, XX.mean %>%
  #                  group_by(SUBJECT) %>% 
  #                  mutate(TP.baseline = min(TP), n = sum(!is.na(value))) %>% 
  #                  filter(TP.baseline == 1) %>% 
  #                  mutate(value.dem = (value - value[TP.baseline])/ mean(value, na.rm = T) * 100) %>% 
  #                  dplyr::select(-c(TP.baseline, n)) %>% ungroup())

  ######################################
  mycoefs = train_res$prepared_data$mycoefs
  
  common = intersect(mypars$subject, XX.init$SUBJECT)
  
  mypars.common = mypars%>% filter(subject %in% common)%>% arrange
  mycoefs.common =  mycoefs[XX.init$SUBJECT %in% common, ]
  #mycoefs.common =  scores[XX.init$SUBJECT %in% common, ]
  
# test for different parnames
  for(parname in parnames){
    
  print(paste(parname, label))
  
  # select trainers
  # compute correlations
  intervention = mypars.common$group == "Intervention"
  cors.intervention = apply(mycoefs.common[intervention, ], 2, function(x) cor.test(x, mypars.common[intervention, parname])$estimate)
  cors.control = apply(mycoefs.common[!intervention, ], 2, function(x) cor.test(x, mypars.common[!intervention, parname])$estimate)
  cors = apply(mycoefs.common, 2, function(x) cor.test(x, mypars.common[, parname])$estimate)
#  plot(sort(cors), type = 'l', main = paste(fulllabel, parname), ylim = c(-1, 1))
#  lines(sort(cors.intervention), type = 'l', col = 'red')
#  lines(sort(cors.control), type = 'l', col = 'blue')
#  abline(0.05, 0, col = 'red')
  hist(cors.intervention, 50, main = paste(fulllabel, parname), xlim = c(-.5, .5))
  #if (min(cors.intervention) < 0.05) print(cors.intervention)
  
  pred = scores %*% as.matrix(colMeans(train_res$coefs.iter[, valid.features]))
  cor.pred = cor.test(pred[XX.init$SUBJECT %in% common][intervention], mypars.common[intervention, parname])
  #if (cor.pred$p.value < 0.05) print(cor.pred)
  
}
}

# X.sem = XX.all %>% group_by(TP, GROUP, FULL_LABEL) %>% dplyr::summarise(y.mean = mean(value.dem, na.rm = T),
#                                                        y.sem = se(value.dem)) %>% ungroup()

# myplot =  ggplot(NULL) +
#   geom_line(data = XX.all, aes(
#     x = TP,
#     group = SUBJECT,
#     col = GROUP, 
#     y = value.dem
#   ), size = 0.1, alpha = 0.6) +
#   geom_line(data = X.sem, aes(
#     x = TP,
#     group = GROUP,
#     col = GROUP,
#     y = y.mean
#   )) +
#   geom_point(data = X.sem, aes(
#     x = TP,
#     group = GROUP,
#     col = GROUP,
#     y = y.mean
#   )) +
#   geom_errorbar(data = X.sem,
#                 aes(
#                   x = TP,
#                   ymax = y.mean + y.sem,
#                   ymin = y.mean - y.sem,
#                   group = GROUP,
#                   col = GROUP
#                 )) +
#   scale_colour_manual(values = myPalette) + theme_lh + 
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.title = element_blank(), legend.position = "bottom") +
#   xlab('Scanning session') + 
#   ylab(MYYLABEL) + 
#   ylim(-10, 10) + 
#   facet_grid(. ~ FULL_LABEL)
# #+ 
# #  scale_x_continuous(breaks = seq(7)) + 
# #  ggtitle(title) + geom_hline(yintercept = 0, size = 0.2)
# 
# print(myplot)

#ggsave(file.path(FIGS_DIR, 'Timecourse.png'), dpi = DPI)

NCOL = 4
ggsave(
  filename = file.path(FIGS_DIR, 'Timecourse.png'),
  plot = ggarrange(
    plotlist = myplots,
    nrow = ceiling(length(myplots) / NCOL),
    ncol = min(NCOL, length(myplots))
  ),
  width = 21,
  height = 13.2
)
