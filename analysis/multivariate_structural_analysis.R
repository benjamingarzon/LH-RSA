rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(pracma)
library(freesurfer)
library(caret)
library(glmnet)
library(forecast)
library(e1071)
library(doParallel)
library(foreach)

set.seed(1)
#sync permutations

# can put back system as covariate
# project to sphere
# pca
# fit models
# null distribution

# make sure variances are ok xx
# project weights xx
# iterate correlations
# make sure that null distribution is correct
# remove 7 timepoint 
# screeplot
# RBM
# remove outliers

# try on whole brain xx
# try repeated cv xx
# output results xx
# balance classes xx
# cross-validation... xx
#tuning xx
# permutation xx
# remove effect of classic... xx
# parallelize xx

cat("\014")

source('~/Software/ImageLMMR/ImageLMMR.R')
setwd("~/Software/LeftHand/analysis")

SUBJECTS_DIR= '~/Data/LeftHand/Lund1/freesurfer/'

OUTPUTDIR='multivariate-tess'
NPROCS = 10

source('./load_covariates.R')

###############################################################
# Define the pipeline
###############################################################

open_data = function(IMAGING_FILE, MASK_FILE){
  
  imaging <- fast_readnii(IMAGING_FILE)
  mask <- fast_readnii(MASK_FILE)
  
  shape <- dim(imaging)
  n_scans <- shape[4]
  n_voxels <- sum(mask > 0)
  print(paste('Voxels in mask:', n_voxels))  
  
  imaging.mat <- matrix(0, n_scans, n_voxels) 
  
  #get_time()
  print('Reshaping data')  
  for (i in seq(n_scans)) {
    #print(i)
    imaging.mat[i, ] <- imaging[, , , i][mask > 0]
  }
  
  rm(imaging)
  
  imagdim = dim(imaging.mat)
  print(paste('Final imaging data size', imagdim[1], imagdim[2]))
  
  #get_time()
  
  return(imaging.mat)
}

prepare_data.old = function(X, DATA, permutate = F){

  XX = cbind(DATA[c("SUBJECT", "TP", "GROUP.NUM")], X)
  XX.melt = melt(XX, id.vars = c("SUBJECT", "TP", "GROUP.NUM"), variable_name = "VERTEX")
  print(colnames(XX.melt))
  browser()
  XX = cast(XX.melt, SUBJECT + GROUP.NUM + VERTEX ~ TP )
  
  complete = rowSums(is.na(XX)) == 0
  XX.init <- XX[, -c(1, 2, 3)]
  incomplete = which(!complete)
  #get_time()
  print("Imputing")
  print(dim(XX.init))
  for (i in incomplete){
    uu = as.numeric(na.interp(ts(t(XX.init[i, ]))))
    if (!permutate) XX.init[i, ] = uu
    else XX.init[i, ] = sample(uu)
  }
  XX[, -c(1, 2, 3)] = XX.init

  print("Reshaping")
  #get_time()
  XX.melt = melt(XX, id.vars = c("SUBJECT", "GROUP.NUM", "VERTEX"), variable_name = "TP")
  XX.cast = cast(XX.melt, SUBJECT + GROUP.NUM ~ TP + VERTEX) %>% arrange(SUBJECT)
  XX.init = XX.cast[, -c(1,2)]
  XX.demeaned = XX.init - rowMeans(XX.init)
  
  y = XX.cast$GROUP.NUM
  #  if (permutate) y = sample(y)
  
  return(list(X = XX.demeaned, y = y))
}


prepare_data = function(X, DATA){

  mycoefs = NULL
  myvars = c("SUBJECT", "TP", "TRAINING", "TRAINING.Q", "GROUP.NUM")
  XX = cbind(DATA[myvars], X)
  XX.melt = reshape::melt(XX, id.vars = myvars, variable_name = "VERTEX")
  XX.init = unique(XX.melt[c("SUBJECT", "GROUP.NUM")])
  print("Fitting lmer model")
  which.model = NULL
  #browser()
  for (subject in XX.init$SUBJECT){
    
      
      mydata = subset(XX.melt, SUBJECT == subject)
      #model.b = brm(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING + TRAINING.Q |VERTEX), data = mydata) 
      
      #                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      model.q = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING + TRAINING.Q |VERTEX), data = mydata) 
      model.q2 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING.Q |VERTEX), data = mydata) 
      model.l = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING |VERTEX), data = mydata)
      model.0 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 |VERTEX), data = mydata)
      #myBIC = BIC(model.0, model.l, model.1)
      mymodels = list(c = model.0, l = model.l, q = model.q, q2 = model.q2)
      is.sing = sapply(mymodels, isSingular)
      mymodels = mymodels[!is.sing]
      myBIC = sapply(mymodels, BIC)
      which.model[subject] = names(which.min(myBIC))
      model = mymodels[[which.model[subject]]]
      RFX = ranef(model)$VERTEX
      FFX = fixef(model)
      #mycoefs = rbind(mycoefs, c(XX.coefs[, 2], XX.coefs[, 3]))
      if (which.model[subject] == 'c') RFXX = cbind(RFX, 0, 0)
      if (which.model[subject] == 'l') RFXX = cbind(RFX, 0)
      if (which.model[subject] == 'q2') RFXX = cbind(RFX[, 1], 0, RFX[ ,2])
      XX.coefs = t(t(RFXX) + FFX) 
      mycoefs = rbind(mycoefs, c(XX.coefs[, 2], XX.coefs[, 3])) # ignore  intercept
  }
  
  global.model = min(which.model)
  if (global.model == 1) mycoefs = NULL 
  if (global.model == 2) mycoefs = mycoefs[ , 1:(ncol(mycoefs)/2)]

  y = XX.init$GROUP.NUM
  return(list(X = mycoefs, y = y, which.model = which.model))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

get_pca = function(X, NREPS = 101){
  print("Doing PCA...")

  sdev = NULL
  for (i in seq(NREPS)){
    # shuffle columns and do PCA
    if (i == 1) XX = X else XX = apply(X, 2, sample, replace = F)
    pca.0 = prcomp(XX, center = F, scale = F)
    sdev = rbind(sdev, apply(pca.0$x, 2, sd))
    if (i == 1) pca = pca.0 # keep this one 
  }
  
  # select how many components to retain
  quants = apply(sdev, 2, function(x) quantile(x[-1], 0.01))
  comps = sdev[1, ] > quants
  ncomps = max(which(comps)) + 1
  plot(sdev[1, ])
  lines(quants, col = "red")
  pca$x = pca$x[, 1:ncomps, drop = F]
  pca$rotation = pca$rotation[, 1:ncomps, drop = F]
  pca$ncomps = ncomps
  print(paste(ncomps, "components"))
  return(pca)
}


train_model = function(X, y, permutate = F, folds){
  TUNELENGTH = 3 #5
  NFOLDSINNER = 3 #10
  NREPEATSINNER = 10
  NITER = 2 # 20 
  
  #mycv = trainControl(method = "repeatedcv", repeats = NREPEATSINNER, number = NFOLDSINNER)
  mycv = trainControl(method = "cv", number = NFOLDSINNER)
  
  # standardize
  #print(dim(X))
  #pca = prcomp(X, center = F, scale = F)
  pca = get_pca(X)
  
  print(dim(pca$x))
  #X.rec = t(t(pca$x %*% t(pca$rotation))) # || t(t(pca$x %*% t(pca$rotation)) * pca$scale + pca$center)
  YY = as.data.frame(pca$x)
  
  y = y + 1 # between 1 and 2
  
  y.pred = y.test = coefs.iter = NULL
  print("Doing classification")
  
  for (fold in folds){
    y.train = as.factor(y[fold])
    print(paste("Fold", fold))
    if (permutate) y.train = sample(y.train)
    
    train_set = YY[fold, ] 
    test_set = YY[-fold, ]
    print(dim(train_set))
    
    pvalues = apply(train_set, 2,  function(x) t.test(x ~ y.train)$p.value)
    myfeatures = which(pvalues < 0.1)
#    if (length(myfeatures) <= 1 ) myfeatures = unique(c(myfeatures, sample(seq(length(pvalues)), size = 2, replace = F))) # get a couple
#    preproc = preProcess(train_set[myfeatures], method = c("center", "scale"))
#    train_set = predict(preproc, train_set[, myfeatures, drop = F])
#    test_set = predict(preproc, test_set[, myfeatures, drop = F])

    preproc = preProcess(train_set, method = c("center", "scale"))
    train_set = predict(preproc, train_set)
    test_set = predict(preproc, test_set)
    
    smax = max(table(y.train))
    
    y.pred.iter = NULL
    
    for (iter in 1:NITER){
      
      # resample and balance classes
      s1 = sample(which(y.train == 1), size = smax, replace = T)
      s2 = sample(which(y.train == 2), size = smax, replace = T)
      s = c(s1, s2)
      
      train.s = train_set[s, ]
      y.s = y.train[s]
#      tryCatch({
        mysvm = train(train.s, y.s,
                      method = "svmLinear", 
                      trControl = mycv, 
                      tuneGrid = expand.grid( C = 10^seq(-2, 0.5, length.out = TUNELENGTH)),
                      verbose = F)
        coefs <- mysvm$finalModel@coef[[1]]
        mat <- mysvm$finalModel@xmatrix[[1]]
        
        #mycoefs = pvalues*0
        #mycoefs[myfeatures] = coefs %*% mat
        mycoefs = coefs %*% mat
        coefs.iter = rbind(coefs.iter, mycoefs)
        y.pred.iter = rbind(y.pred.iter, predict(mysvm, newdata = test_set))
        print(paste(iter, mysvm$bestTune), sep = ':')
        
      # }, error = function(error_message) { 
      #   message("We have a problem training")
      #   message(error_message) 
      #   save(train.s, y.s, file = "/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results/tests/multivariate-tess/test.rds")
      #   return(NA)
      # }
      # ) # try 
      # 
      #plot(mysvm$results$C, mysvm$results$Accuracy, type = 'l')
    }
    y.pred = c(y.pred, apply(y.pred.iter, 2, Mode))
    #print(y.pred)
    y.test = c(y.test, y[-fold])
  }
  cm = confusionMatrix(data = as.factor(y.pred), reference = as.factor(y.test))
  accuracy = cm$byClass['Balanced Accuracy']
  print("Finished prediction")
  return(list(accuracy = accuracy, coefs.iter = coefs.iter, pca = pca, y = y))
  
}

analysis = function(label, imaging.mat, DATA, annotation, folds, y, permutate = F){
  print("Start analysis")
  print(label)
  w = which(annotation$colortable$label == label)
  code = annotation$colortable$code[w]  
  mask = annotation$label == code

  X = imaging.mat[, mask]
  
#    if (min(LL$which.model == 2)) browser()
  LL = prepare_data(X, DATA)
  
  train_res = train_model(LL$X, LL$y, permutate, folds)
  if (!exists(train_res)) train_res = list(success = F)
  else train_res$success = T
  
  #print(y == LL$y)
  
  # project data to original space
  valid.features = which(apply(train_res$coefs.iter, 2, function(x) max(table(sign(x))[c('-1','1')], na.rm = T)/length(x) > 0.8))
  if (length(valid.features) > 0){
    scores = train_res$pca$x[, valid.features, drop = F]
    loadings = train_res$pca$rotation[, valid.features, drop = F]
    c1 = scores[train_res$y == 1, ]
    c2 = scores[train_res$y == 2, ]
    x1 = loadings %*% t(c1)
    x2 = loadings %*% t(c2)
    m1 = rowMeans(x1)
    m2 = rowMeans(x2)
    s1 = apply(x1, 1, sd)
    s2 = apply(x2, 1, sd)
    
    #unfold in time
    
    indices = t(sapply(rownames(loadings), function(x) unlist(strsplit(x, '_'))))
    colnames(indices) = c("TP", "VERTEX")
    data.proj = as.data.frame(cbind(indices, m1 = m1, m2 = m2, s1 = s1, s2 = s2)) %>% mutate_all(as.character) %>% mutate_all(as.numeric)
    train_res$c1 = c1
    train_res$c2 = c2
    train_res$s1 = s1
    train_res$s2 = s2
    train_res$data.proj = data.proj
    train_res$scores = scores
    train_res$loadings = loadings
  }
  train_res$mask = mask
  train_res$valid.features = valid.features
  train_res$label = label
  return(list(train_res = train_res))
}


doit = function(WD, 
                MASK = 'rh.cortex.mask.nii.gz', 
                IMAGES_LIST = 'rh.thickness.txt',
                IMAGING_NAME = 'rh.thickness.nii.gz',
                ANNOT_FILE = NULL,
                LABEL_FILE = NULL,
                ACC_FILE = NULL,
                NPERMS = 0
                #                to_gifti = '', 
                #                NPERMS = 0, 
                #                alpha = 0.05,
                #                shuffle_by = NULL, remove_outliers = T
)
{
  
  setwd(WD)
  IMAGES = read.table(file.path(WD, IMAGES_LIST))
  IMAGING_FILE = file.path(WD, IMAGING_NAME)
  MASK_FILE = file.path(WD, MASK)
  
  if (is.null(IMAGES$V3)) {
    DATA = data.frame(
      SUBJECT = IMAGES$V1, 
      TP = IMAGES$V2
    )
    DATA = merge(DATA, covars.table, by = c("SUBJECT", "TP"), all.x = T) %>% arrange(SUBJECT, TP)
    
  } 
  else {
    DATA = data.frame(
      SUBJECT = IMAGES$V1, 
      TP = IMAGES$V2,
      DEPTH = 1 - IMAGES$V3 # depth starting from cortical surface/ inverse depth
    )
    DATA = merge(DATA, covars.table, by = c("SUBJECT", "TP"), all.x = T) %>% arrange(-DEPTH, SUBJECT, TP)
    
  }
  
  View(DATA)
  print(DATA[!complete.cases(DATA), ])
  
  imaging.mat = open_data(IMAGING_FILE, MASK_FILE)
  # exclude also MTX8 
  imaging.mat = imaging.mat[DATA$SYSTEM == "Classic", ]
  DATA = DATA[DATA$SYSTEM == "Classic", ]
  
  excluded = which((DATA$SUBJECT %in% names(which(table(DATA$SUBJECT) < 4))) | !complete.cases(DATA)) 
  
  # remove bad cases 
  if (!is.null(excluded)){
    imaging.mat = imaging.mat[-excluded, ]
    DATA = DATA[-excluded, ]
  }
  
  annotation = read_annotation(ANNOT_FILE)
  labels = unlist(read.table(LABEL_FILE))[-c(1, 2)]
  
  # if (F) {
  #   accuracy = NULL
  #   for (label in labels) { 
  #     print(label)
  #     w = which(annotation$colortable$label == label)
  #     code = annotation$colortable$code[w]
  #     mask = annotation$label == code
  #     
  #     X = imaging.mat[, mask]
  #     LL = prepare_data(X, DATA)
  #     accuracy[label] = train_model(LL$X, LL$y)
  #     print(accuracy)
  #   }
  # } else {
  #   
    cl <- makeCluster(NPROCS)
    registerDoParallel(cl)
    accuracies = NULL
    NFOLDS = 5 # 20
    NREPEATS = 10  # 10
    y = unique(DATA[c("SUBJECT", "GROUP.NUM")]%>%arrange(SUBJECT))$GROUP.NUM
    
    folds = createMultiFolds(y, k = NFOLDS, times = NREPEATS) 
    #output = NULL
    myperms = seq(NPERMS + 1)
    iterations = expand.grid(perm = myperms, label = labels) # reorder to go for first
    
    k = 1
    niters = nrow(iterations)
    accuracies = matrix(0, length(labels), length(myperms))
    colnames(accuracies) = c("value", paste0("perm_", sprintf("%03d", seq(NPERMS))) )

    results = list()
    
    while (T){
      indices = seq(NPROCS*(k - 1) + 1, min(niters, NPROCS*k))
      print(indices)
      if (indices[1] > niters) break
      
      get_time()
      train_res_list = foreach(index = indices, # ###label = labels,
                               .packages = c('caret', 'e1071', 'forecast', 'reshape', 'dplyr', 'lme4', 'brms'),
                               .export = c('prepare_data', 'open_data', 'train_model', 'get_time', 'analysis', 'Mode', 'get_pca'),
                               .combine = 'c') %do% analysis(iterations$label[index], imaging.mat, DATA, annotation, folds, y, 
                                                             permutate = iterations$perm[index] > 1)

      #get accuracies and print out
      print("Finished indices")
      accuracy_list = unlist(lapply(train_res_list, function(x) x$accuracy))
      
      for (index in indices){
        accuracies[iterations$label[i], iterations$perm[i]] = accuracy_list[i]
      }
      
      output = accuracies 
      output$labels = labels
      write.table(output, file = ACC_FILE, sep = ',', col.names = T, row.names = F)
      
      # no need to generate every time but quick
      # binomial distribution
      binomsize = NREPEATS*length((y))
      
      output = matrix(rbinom(length(labels)*NPERMS, binomsize, p = 0.5)/binomsize, length(labels), NPERMS)
      colnames(output) = c("value", paste0("perm_", sprintf("%03d", seq(NPERMS))) )
      output$labels = labels
      write.table(output, file = paste0(ACC_FILE, '.binom'), sep = ',', col.names = T, row.names = F)
      
      results$train_res_list = c(results$train_res_list, train_res_list)
      
      k = k + 1
      
    } # while
    
    # get accuracies
    
    
      #print(accuracies)
#      if(myperm > 1){
#        output = as.data.frame(cbind(as.character(labels), accuracies))
#        colnames(output) = c("label", "value", paste0("perm_", sprintf("%03d", seq(myperm - 1))) )
#      } else {
#        output = as.data.frame(cbind(as.character(labels), accuracies))
#        colnames(output) = c("label", "value")
#      }
      #write.table(output, file = ACC_FILE, sep = ',', col.names = T, row.names = F)
#      get_time()
#    }
#    write.table(output, file = ACC_FILE, sep = ',', col.names = T, row.names = F)
    stopCluster(cl)
    
  results$accuracies = accuracies
  results$data = DATA
  #  save(results, file = file.path(OUTPUT_DIR, 'results.rda'))
  return(results)
}

###############################################################
# Do it
###############################################################

rh.cortex.mask = 'rh.cortex.mask.nii.gz'
lh.cortex.mask = 'lh.cortex.mask.nii.gz'

shuffle_by = c('BETWEEN', 'WITHIN')

alphavoxel = 0.05
alphavertex = 0.025


#####################
# thickness
#####################

DATADIR='/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'
#annotation.rh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations/GlasserParc/rh.HCP-MMP1.annot'
#labels.rh = '/home/benjamin.garzon/Software/LeftHand/masks/rh.parcels.txt'
#annotation.lh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations/GlasserParc/lh.HCP-MMP1.annot'
#labels.lh = '/home/benjamin.garzon/Software/LeftHand/masks/lh.parcels.txt'

annotation.rh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations/fs_LR_32/Icosahedron-162.fsaverage.R.annot'
labels.rh = '/home/benjamin.garzon/Software/LeftHand/masks/rh.tessellation162.txt'
annotation.lh = '/home/benjamin.garzon/Data/LeftHand/Lund1/parcellations//fs_LR_32/Icosahedron-162.fsaverage.L.annot'
labels.lh = '/home/benjamin.garzon/Software/LeftHand/masks/lh.tessellation162.txt'

NPERMS = 5
results.multivariate.rh = doit(DATADIR,
                               MASK = rh.cortex.mask,
                               IMAGES_LIST = 'rh.thickness.txt',
                               IMAGING_NAME = 'rh.thickness.nii.gz',
                               ANNOT_FILE = annotation.rh,
                               LABEL_FILE = labels.rh,
                               ACC_FILE = file.path(DATADIR, 'tests', OUTPUTDIR, 'rh.accuracy'),
                               NPERMS = NPERMS) # remove outliers

results.multivariate.lh = doit(DATADIR,
                               MASK = lh.cortex.mask,
                               IMAGES_LIST = 'lh.thickness.txt',
                               IMAGING_NAME = 'lh.thickness.nii.gz',
                               ANNOT_FILE = annotation.lh,
                               LABEL_FILE = labels.lh,
                               ACC_FILE = file.path(DATADIR, 'tests', OUTPUTDIR, 'lh.accuracy'),
                               NPERMS = NPERMS) # remove outliers

save(results.multivariate.lh, 
     results.multivariate.rh, 
     file = file.path(DATADIR, 'tests', OUTPUTDIR, 'results.rds'))
# collect and plot results

accuracy.rh = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'rh.accuracy'), header = T, sep = ',')
accuracy.lh = read.table(file.path(DATADIR, 'tests', OUTPUTDIR, 'lh.accuracy'), header = T, sep = ',')
inflated.rh = freesurfer_read_surf('/usr/local/freesurfer/subjects/fsaverage/surf/rh.inflated')
inflated.lh = freesurfer_read_surf('/usr/local/freesurfer/subjects/fsaverage/surf/lh.inflated')
accuracy.rh$hemi = "rh"
accuracy.lh$hemi = "lh"

accuracy = rbind(accuracy.rh, accuracy.lh)
means = colMeans(dplyr::select(accuracy, -c(label, hemi)))
maxs = apply(dplyr::select(accuracy, -c(label, hemi)), 2, max)

pvals = apply(dplyr::select(accuracy, -label), 1, function(x) mean(x[1] <= x))
pvals.corr = sapply(accuracy$value, function (x) mean(x <= maxs)) 
res = cbind(dplyr::select(accuracy, c(label, hemi, value)), pvals = pvals, pvals.corr = pvals.corr )%>% arrange(pvals.corr)
print(means)
print(maxs)
res.significant = subset(res, pvals < 0.05)
res.corr = subset(res, pvals.corr < 0.05)

data.max = expand.grid(label = accuracy$label, m = maxs[-1]) %>% filter(label %in% res.significant$label)
data.max = merge(data.max, res.significant, by = 'label')
myplot = ggplot() + 
  geom_violin(data = data.max, aes (x = reorder(label, -value), y = m)) +  
  geom_point(data = res.significant, aes(x = label, y = value, col = hemi), size = 2) + 
  ylim(0.5, 0.8) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+ facet_grid(. ~ hemi)
print(myplot)

write.table(filter(res.corr, hemi == "rh") %>% dplyr::select(c(label, value)), file.path(DATADIR, 'tests', OUTPUTDIR, 'rh.accuracy.results.csv'), 
            col.names = T, sep = ',', row.names = F)
write.table(filter(res.corr, hemi == "lh") %>% dplyr::select(c(label, value)), file.path(DATADIR, 'tests', OUTPUTDIR, 'lh.accuracy.results.csv'), 
            col.names = T, sep = ',', row.names = F)
write.table(dplyr::select(res.corr, c(label, value)), file.path(DATADIR, 'tests', OUTPUTDIR, 'accuracy.csv'), 
            col.names = T, sep = ',', row.names = F)

train_res_list = c(results.multivariate.rh$train_res_list, results.multivariate.lh$train_res_list)
data.proj.all = NULL
for (label in res.corr$label){
  print(label)
  w = which(sapply(train_res_list, function(x) x$label == label))
  train_res = train_res_list[[w]]
  mask.indices = which(train_res$mask)
  data.proj = train_res$data.proj
  data.proj$INDEX = mask.indices[data.proj$VERTEX]
  coords = inflated.rh$vertices
  data.proj$x = coords[data.proj$INDEX, 1] 
  data.proj$y = coords[data.proj$INDEX, 2]
  data.proj$z = coords[data.proj$INDEX, 3]
  data.proj$LABEL = label
  data.proj.all = rbind(data.proj.all, data.proj)
}

# visualize patches 
data.proj.melt = melt(data.proj.all, id.vars = c("TP", "VERTEX", "INDEX", "x", "y", "z", "LABEL"), measure.vars = c("m1", "m2"), variable_name = "GROUP")
myplot = ggplot(data.proj.melt, aes ( x = y, y = z, col = value)) + geom_point() + scale_colour_viridis_c() + facet_grid(LABEL + GROUP ~ TP)
myplot = ggplot(data.proj.melt, aes ( x = TP, y = value, col = GROUP, group = VERTEX)) + geom_line(size = 0.1) + facet_grid(GROUP ~ LABEL)
print(myplot)


BINS = seq(0.2, 0.8, length.out = 30)
hist(accuracy$value, breaks = BINS, xlim = c(0.2, 0.8), ylim = c(0, 15))
hist(accuracy$perm_002, breaks = BINS, add = T, col = rgb(1, 0, 0, 0.2))
hist(xx, breaks = BINS, add = T, col = rgb(0, 0, 1, 0.2))

hist(rowMeans(dplyr::select(accuracy, -c(label, value))), breaks = BINS, add = T, col = rgb(1, 0, 0, 0.2))

hist(as.numeric(subset(accuracy, label == "R_IP2_ROI")[-1]), breaks = BINS, xlim = c(0.2, 0.8), ylim = c(0, 15))
