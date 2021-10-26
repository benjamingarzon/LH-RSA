# compose figs for structural analyses

rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggpubr)
library(png)
library(lme4)
library(lmerTest)
library(reshape2)
library(reshape)
library(doParallel)
library(foreach)



FIGS_DIR = '~/Data/LeftHand/Lund1/figs/structure'
DPI = 1000
FONT.SIZE = 14
NPROCS = 35
NSAMPLES = 100
NVOXELS = 1000
N = 30 # subjects per group
alpha = 0.001
maxvals = seq(0, 10, 1) / 100

GMV.png = file.path(FIGS_DIR, 'PowerGMV.png')
thickness.png = file.path(FIGS_DIR, 'PowerThickness.png')
doGMV = file.exists(thickness.png)


if (doGMV) {
  DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/cat12'
  output_file = 'PowerGMV.png'
  exclude_legend = F
} else {
  DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'
  load(file.path(DATADIR, 'tests', 'thickness', 'reliability-mask.rh', 'results.rda'))
  
  DATADIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/cat12'
  output_file = 'PowerThickness.png'
  exclude_legend = T
}

load(file.path(DATADIR, 'tests', 'reliability-mask', 'results.rda'))
imaging = results$imaging.mat
data = results$data[-results$excluded,]

indices = sort(sample(seq(ncol(imaging)), NVOXELS, replace = F))
imaging = imaging[data$GROUP == 'Control' &
                    data$SYSTEM == 'Classic', indices]
feature_names = paste("voxel", seq(ncol(imaging)), sep = "")
colnames(imaging) = feature_names
data = subset(data, GROUP == 'Control' & SYSTEM == 'Classic')
data.complete = cbind(data, imaging)

data.sum = data.complete %>% group_by(SUBJECT) %>% summarise_at(feature_names, .funs = list(mean = mean, sd = sd))
#sample subject, compute mean and sd vector and generate data

# define some functions

mymodel = function(x, y, maxval, type = 'linear') {
  n = length(x) - 1
  xx = (x - 1) / n
  xx2 = (xx - mean(xx)) ^ 2
  xx2 = 1 - xx2 / max(xx2)
  xxa = xx2
  xxa [(n / 2 + 1):(n + 1)] = 1
  if (type == 'Linear')
    z = xx * maxval * mean(y) + y
  if (type == 'Asymptotic')
    z = xxa * maxval * mean(y) + y
  if (type == 'Quadratic')
    z = xx2 * maxval * mean(y) + y
  return(z)
}

sample_subjects = function(data.sum, N, NTPS, maxval, model_type) {
  means = data.sum[, grep ('_mean', colnames(data.sum))]
  sds = data.sum[, grep ('_sd', colnames(data.sum))]
  
  subjects = data.sum$SUBJECT
  sampled_subjects = sample(subjects, N * 2, replace = T)
  
  subjects = unique(data$SUBJECT)
  data.new = NULL
  for (subject in sampled_subjects) {
    data.new = rbind(data.new,
                     mapply(function(x, y)
                       rnorm(NTPS, x, y),
                       means[subjects == subject,],
                       sds[subjects == subject,]))
  }
  
  colnames(data.new) = feature_names
  data.new = as.data.frame(data.new)
  data.new$SUBJECT = as.numeric(sapply(seq(2 * N), function(x)
    rep(x, NTPS)))
  data.new$GROUP = c(rep("Control", nrow(data.new) / 2), rep("Intervention", nrow(data.new) /
                                                               2))
  data.new$TP = rep(seq(NTPS), 2 * N)
  for (var in feature_names) {
    dd = data.new %>% filter(GROUP == 'Intervention') %>% group_by(SUBJECT) %>%
      mutate(value = mymodel(TP, get(var), maxval, type = model_type))
    dd[var] = dd$value
    dd$value = NULL
    data.new[data.new$GROUP == 'Intervention',] = dd
  }
  
  #plot.sim = ggplot(data.new, aes(x = TP, y = voxel20, group = SUBJECT)) + geom_line()
  #plot.true = ggplot(data.complete, aes(x = TP, y = voxel20, group = SUBJECT)) + geom_line()
  #print(ggarrange(plot.true, plot.sim))
  
  return(data.new)
}

fit_model = function(var, data.sim, model_type){

  if (model_type == 'Linear'){ 
    myformula = paste(var, '1 + GROUP*TP.1 + (1|SUBJECT)', sep =  '~')
    model = lmer(as.formula(myformula), data = data.sim)
    pval = summary(model)$coefficients['GROUPIntervention:TP.1', 'Pr(>|t|)']
  }
  
  if (model_type == 'Asymptotic') {
    myformula = paste(var, '1 + GROUP*TP.a + (1|SUBJECT)', sep =  '~')
    model = lmer(as.formula(myformula), data = data.sim)
    pval = summary(model)$coefficients['GROUPIntervention:TP.a', 'Pr(>|t|)']
  }
  
  if (model_type == 'Quadratic') {
    myformula = paste(var, '1 + GROUP*TP.2 + (1|SUBJECT)', sep =  '~')
    model = lmer(as.formula(myformula), data = data.sim)
    pval = summary(model)$coefficients['GROUPIntervention:TP.2', 'Pr(>|t|)']
  }
  
  
  return(pval)

}

get_power = function(data.sum, N, NTPS, maxval, model_type = 'linear') {
  pvals = NULL
  cl <- makeCluster(NPROCS)
  registerDoParallel(cl)
  
  for (k in seq(NSAMPLES)) {
#    print(k)
    data.sim = sample_subjects(data.sum, N, NTPS, maxval, model_type)
    m = max(data.sim$TP)
    data.sim$TP.1 = (data.sim$TP - 1)/(m - 1)
    data.sim$TP.2 = -(data.sim$TP.1 - 0.5)^2
    data.sim$TP.2 = (data.sim$TP.2 - min(data.sim$TP.2))*4
    data.sim$TP.a = data.sim$TP.2
    data.sim$TP.a[data.sim$TP> (m-1)/2] = 1

    pval = rep(0, length(feature_names))
    names(pval) = feature_names

    pval = foreach(var = feature_names,
                             .packages = c('lmerTest','lme4'),
                             .export = 'fit_model',
                             .combine = 'c') %dopar% fit_model(var, data.sim, model_type)

    pvals = rbind(pval, pvals)
  }
  
  stopCluster(cl)
  power = colMeans(pvals < alpha)
  return(mean(power))
}

proportions = c(
  astrocytes = 9,
  collaterals = 29,
  dendrites = 26,
  soma = 11,
  oligodendrocytes = 2.5,
  synapses = 5.7,
  total = 100
) / 100

simulate = function(NTPS, model_type = 'linear') {
  powers = NULL
  for (maxval in maxvals) {
    print(maxval)
    powers = c(powers, get_power(data.sum, N, NTPS, maxval, model_type))
  }
  mymat = as.data.frame(sapply(proportions, function(x)
    maxvals / x)) * 100
  mymat$power = powers
  mymat = reshape2::melt(
    mymat,
    id.vars = 'power',
    variable.name = 'type',
    value.name = 'fraction'
  )
  mymat$model_type = model_type
  mymat$NTPS = NTPS
  return(mymat)  
}




mymat.linear.7 = simulate(7, 'Linear')
mymat.asymptotic.7 = simulate(7, 'Asymptotic')
mymat.quadratic.7 = simulate(7, 'Quadratic')
mymat.linear.3 = simulate(3, 'Linear')
mymat.asymptotic.3 = simulate(3, 'Asymptotic')
mymat.quadratic.3 = simulate(3, 'Quadratic')

mymat = rbind(mymat.linear.7, 
              mymat.asymptotic.7,
              mymat.quadratic.7,
              mymat.linear.3,
              mymat.asymptotic.3,
              mymat.quadratic.3)
mymat$type_NTPS= paste(mymat$type, mymat$NTPS)
plot.power = ggplot(mymat, aes(
  x = fraction,
  y = power*100,
  col = type,
  group = type_NTPS,
  linetype = as.factor(NTPS)
)) +
  geom_line() + geom_point() +
  coord_cartesian(xlim = c(0, 30)) + theme_classic() + 
  xlab('Necessary increase (%)') + ylab('Power (%)') + labs(col = 'Cellular constituent', linetype = 'Number of sessions') +
  geom_hline(yintercept = 80, linetype = 2) + facet_grid(. ~ model_type) + theme_classic() + theme(strip.background = element_blank(), 
                                                                                                    legend.position = 'top', legend.box = 'vertical',
                                                                                                    text = element_text(size = FONT.SIZE)) 
if (exclude_legend) plot.power = plot.power + theme(legend.position = 'none')
print(plot.power)

ggsave(filename = file.path(FIGS_DIR, output_file), plot = plot.power, dpi = DPI, width = 16, units = 'cm') 

if (file.exists(thickness.png) & file.exists(GMV.png)){
img.GMV = readPNG(GMV.png)
img.GMV = ggplot() + background_image(img.GMV)
img.thickness = readPNG(thickness.png)
img.thickness = ggplot() + background_image(img.thickness)

plot.final = ggarrange(img.thickness, img.GMV, 
                             labels = c("A", "B"), 
                             ncol = 1, nrow = 2, font.label = list(size = FONT.SIZE))

ggsave(file.path(FIGS_DIR, 'Power.png'), plot = plot.final, dpi = DPI, width = 16, height = 26)

}
