###############################################################
# Structural analysis pipelines
###############################################################

se = function(x) sd(x, na.rm = T)/sqrt(length(x))

###############################################################
# Univariate
###############################################################
structural_analysis_univariate = function(WD, MYTEST, OD, 
                MASK = 'rh.thickness.mask.nii.gz', 
                IMAGES_LIST = 'rh.thickness.txt',
                IMAGING_NAME = 'rh.thickness.nii.gz',
                to_gifti = '', 
                NPERMS = 0, 
                alpha = 0.05,
                shuffle_by = NULL, upsample = NULL, remove_outliers = F)
{
  
  setwd(WD)
  IMAGES = read.table(file.path(WD, IMAGES_LIST))
  IMAGING_FILE = file.path(WD, IMAGING_NAME)
  MASK_FILE = file.path(WD, MASK)
  
  OUTPUT_DIR = file.path(WD, OD)
  
  # create matrix of regressors
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
  
  # exclude subject if less than 4 timepoints
  excluded = which((DATA$SUBJECT %in% names(which(table(DATA$SUBJECT) < 4))) | !complete.cases(DATA))
  if ( NPERMS == 0)
    results = vbanalysis(
      IMAGING_FILE,
      OUTPUT_DIR, 
      DATA, 
      MASK_FILE,
      MYTEST,
      remove_outliers = remove_outliers, 
      excluded = excluded, 
      to_gifti = to_gifti, alpha = alpha, upsample = upsample
    )
  else 
    results = vbanalysis_perm(
      IMAGING_FILE,
      OUTPUT_DIR, 
      DATA, 
      MASK_FILE,
      MYTEST,
      remove_outliers = remove_outliers, 
      excluded = excluded, 
      to_gifti = to_gifti, alpha = alpha, upsample = upsample,
      NPERMS = NPERMS,
      shuffle_by = shuffle_by
    )
  
  results$data = DATA
  save(results, file = file.path(OUTPUT_DIR, 'results.rda'))
  return(results)
}

###############################################################
# Multivariate
###############################################################

open_data = function(IMAGING_FILE, MASK_FILE) {
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
    imaging.mat[i,] <- imaging[, , , i][mask > 0]
  }
  
  rm(imaging)
  
  imagdim = dim(imaging.mat)
  print(paste('Final imaging data size', imagdim[1], imagdim[2]))
  
  #get_time()
  
  return(imaging.mat)
}

prepare_data = function(X, DATA) {
  mycoefs = NULL
  myvars = c("SUBJECT",
             "TP",
             "FD", 
             "SYSTEM",
             "TRAINING",
             "TRAINING.Q",
             "TRAINING.C",
             "GROUP.NUM")
  XX = cbind(DATA[myvars], X)
  XX.melt = reshape::melt(XX, id.vars = myvars, variable_name = "VERTEX")
  XX.init = unique(XX.melt[c("SUBJECT", "GROUP.NUM")])
  
  # regress out motion and SYSTEM
  
  mymod = lm(value ~ 1 + FD + SYSTEM + TRAINING + TRAINING.Q, data = XX.melt)
  XX.melt$value = predict(mymod, XX.melt %>% mutate(TRAINING = 0, TRAINING.Q = 0, FD = 0, SYSTEM = 'Classic')) + resid(mymod)
  
  print("Fitting lmer model")
  which.model = NULL
  for (subject in XX.init$SUBJECT) {
    
    mydata = subset(XX.melt, SUBJECT == subject)

    if (length(unique(mydata$SYSTEM)) > 1){
    model.q2 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING + TRAINING.Q |
                                                           VERTEX),
                    data = mydata)
    model.q1 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING.Q |
                                                           VERTEX),
                    data = mydata)
    model.l1 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING |
                                                           VERTEX), data = mydata)
    model.l0 = lmer(value ~ 1 + TRAINING + (1 + TRAINING | VERTEX), data = mydata)
    model.c0 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 | VERTEX), data = mydata)
    
    } else {
    model.q2 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING + TRAINING.Q |
                                                                    VERTEX),
                    data = mydata)
    model.q1 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING.Q |
                                                                    VERTEX),
                    data = mydata)
    model.l1 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 + TRAINING |
                                                                    VERTEX), data = mydata)
    model.l0 = lmer(value ~ 1 + TRAINING + (1 + TRAINING | VERTEX), data = mydata)
    model.c0 = lmer(value ~ 1 + TRAINING + TRAINING.Q + (1 | VERTEX), data = mydata)
    }
        #results.l1 = psycho::analyze(model.l1, CI = 95)
    #, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    
    mymodels = list(
      c0 = model.c0,
      l0 = model.l0,
      l1 = model.l1,
      q1 = model.q1,
      q2 = model.q2
    )
    is.sing = sapply(mymodels, isSingular)
    is.sing[1] = F # if first one is singular, still take it
    mymodels = mymodels[!is.sing]
    myBIC = sapply(mymodels, BIC)
    which.model[subject] = names(which.min(myBIC))
    #        if (which.model[subject] == "c0") browser()
    #        else
    print(paste("Winning model:", which.model[subject]))
    model = mymodels[[which.model[subject]]]

    RFX = ranef(model)$VERTEX
    FFX = fixef(model)
    if (which.model[subject] == 'c0')
      RFXX = cbind(RFX, 0, 0)
    if (which.model[subject] == 'l1')
      RFXX = cbind(RFX, 0)
    if (which.model[subject] == 'l0') {
      RFXX = cbind(RFX, 0)
      FFX = c(FFX, 0)
    }
    if (which.model[subject] == 'q2')
      RFXX = cbind(RFX[, 1:3])
    if (which.model[subject] == 'q1')
      RFXX = cbind(RFX[, 1], 0, RFX[, 2])
    
    XX.coefs = t(t(RFXX) + FFX)
    mycoefs = rbind(mycoefs, c(XX.coefs[, 2], XX.coefs[, 3])) # ignore  intercept
    
  }

  ncoefs = ncol(mycoefs) / 2
  colnames(mycoefs) = unlist(lapply(seq(2), paste, seq(ncoefs), sep = '_'))
  
  y = XX.init$GROUP.NUM
  return(list(
    mycoefs = mycoefs,
    y = y,
    X = X,
    subject = XX.init$SUBJECT,
    DATA = DATA,
    XX.melt = XX.melt,
    which.model = which.model
  ))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

get_pca = function(X, NREPS = 101) {
  print("Doing PCA...")
  
  sdev = NULL
  for (i in seq(NREPS)) {
    # shuffle columns and do PCA
    if (i == 1)
      XX = X
    else
      XX = apply(X, 2, sample, replace = F)
    pca.0 = prcomp(XX, center = F, scale = F)
    sdev = rbind(sdev, apply(pca.0$x, 2, sd))
    if (i == 1)
      pca = pca.0 # keep this one
  }
  
  # select how many components to retain
  quants = apply(sdev, 2, function(x)
    quantile(x[-1], 0.01))
  comps = sdev[1,] > quants
  ncomps = max(which(comps)) + 1
  plot(sdev[1,])
  lines(quants, col = "red")
  pca$x = pca$x[, 1:ncomps, drop = F]
  pca$rotation = pca$rotation[, 1:ncomps, drop = F]
  pca$ncomps = ncomps
  print(paste(ncomps, "components"))
  return(pca)
}


train_model = function(X,
                       y,
                       permutate = 0,
                       folds = NULL,
                       dopca = T) {
  TUNELENGTH = 5
  myTuneGrid = expand.grid(C = 1)#, #expand.grid( C = 10^seq(-2, 0.5, length.out = TUNELENGTH))
  NFOLDSINNER = 10
  NREPEATSINNER = 10
  NITER = 10
  
  #mycv = trainControl(method = "repeatedcv", repeats = NREPEATSINNER, number = NFOLDSINNER)
  mycv = trainControl(method = "cv", number = NFOLDSINNER)
  
  if (dopca) {
    pca = get_pca(X)
    YY = as.data.frame(pca$x)
  }
  else {
    pca = NULL
    YY = X
  }
  y = y + 1 # between 1 and 2
  
  y.pred = y.test = coefs.iter = NULL
  print("Doing classification")
  
  for (fold in folds) {
    y.train = as.factor(y[fold])
    # print(paste("Fold", fold))
    if (permutate > 0) {
      set.seed(permutate)
      y.train = sample(y.train)
    }
    
    train_set = YY[fold,]
    test_set = YY[-fold,]
    
    pvalues = apply(train_set, 2,  function(x)
      t.test(x ~ y.train)$p.value)
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
    if (NITER > 0) {
      for (iter in 1:NITER) {
        # resample and balance classes
        s1 = sample(which(y.train == 1),
                    size = smax,
                    replace = T)
        s2 = sample(which(y.train == 2),
                    size = smax,
                    replace = T)
        s = c(s1, s2)
        
        train.s = train_set[s,]
        y.s = y.train[s]
        mysvm = train(
          train.s,
          y.s,
          method = "svmLinear",
          #trControl = mycv,
          tuneGrid = myTuneGrid,
          verbose = F
        )
        coefs <- mysvm$finalModel@coef[[1]]
        mat <- mysvm$finalModel@xmatrix[[1]]
        
        #mycoefs = pvalues*0
        #mycoefs[myfeatures] = coefs %*% mat
        mycoefs = coefs %*% mat
        coefs.iter = rbind(coefs.iter, mycoefs)
        y.pred.iter = rbind(y.pred.iter, predict(mysvm, newdata = test_set))

        #plot(mysvm$results$C, mysvm$results$Accuracy, type = 'l')
      }
      y.pred = c(y.pred, apply(y.pred.iter, 2, Mode))
    } else {
      # no bagging
      mysvm = train(
        train_set,
        y.train,
        method = "svmLinear",
        #trControl = mycv,
        tuneGrid = myTuneGrid,
        verbose = T
      )
      coefs <- mysvm$finalModel@coef[[1]]
      mat <- mysvm$finalModel@xmatrix[[1]]
      
      mycoefs = coefs %*% mat
      coefs.iter = rbind(coefs.iter, mycoefs)
      y.pred = c(y.pred, predict(mysvm, newdata = test_set))
    }
    #print(y.pred)
    y.test = c(y.test, y[-fold])
  } # across folds
  cm = confusionMatrix(data = as.factor(y.pred), reference = as.factor(y.test))
  accuracy = cm$byClass['Balanced Accuracy']
  print("Finished prediction")
  return(list(
    accuracy = accuracy,
    coefs.iter = coefs.iter,
    mycoefs = mycoefs, 
    pca = pca,
    y = y
  ))
  
}

analysis = function(label,
                    imaging.mat,
                    DATA,
                    annotation,
                    folds,
                    y,
                    hemi, 
                    permutate = 0,
                    procdir = NULL) {
  print("Start analysis")
  print(label)
  w = which(annotation$label == as.character(label))
  mask = annotation$index == annotation$code[w]

  X = imaging.mat[, mask]
  myfile = ifelse(is.null(hemi), 
                  file.path(procdir, paste0(label, '.rds')), 
                  file.path(procdir, paste0(hemi, '.', label, '.rds')))
  if (file.exists(myfile))
    load(myfile)
  else
  {
    prepared_data = prepare_data(X, DATA)
    save(prepared_data, file = myfile)
  }
  train_res = train_model(prepared_data$mycoefs, prepared_data$y, permutate, folds)
  
  if (!exists('train_res'))
    train_res = list(success = F)
  else
    train_res$success = T
  
  #print(y == prepared_data$y)
  if (permutate == 0) train_res$prepared_data = prepared_data
  
  train_res$mask = mask
  train_res$label = label
  train_res$perm = permutate
  if (!is.null(hemi)) train_res$hemi = hemi
  return(list(train_res = train_res))
}


load_annotation = function(ANNOT_FILE, MASK_FILE){
  # surface parcellation
  if (length(ANNOT_FILE) == 1) {
    myannot = read_annotation(ANNOT_FILE)
    mask = fast_readnii(MASK_FILE)
    annotation = list(code = myannot$colortable$code,
                      label = as.character(myannot$colortable$label[mask>0]),
                      index = myannot$label)
  
  } else {
    annotvol = fast_readnii(ANNOT_FILE[1])
    mask = fast_readnii(MASK_FILE)
    lut = read.table(ANNOT_FILE[2])
    colnames(lut) = c("code", "label", "R", "G", "B", "alpha")
    labels = sapply(as.character(lut$label), function(x) gsub("-0*", "-", x))
    index = annotvol[mask>0]
    annotation = list(code = lut$code, label = labels, index = index) 
  }
}

structural_analysis_multivariate = function(WD,
                MASK = 'rh.cortex.mask.nii.gz',
                IMAGES_LIST = 'rh.thickness.txt',
                IMAGING_NAME = 'rh.thickness.nii.gz',
                HEMI = 'rh', 
                ANNOT_FILE = NULL,
                LABEL_FILE = NULL,
                ACC_FILE = NULL,
                NPERMS = 0,
                OUTDIR = NULL)
{
  setwd(WD)
  IMAGES = read.table(file.path(WD, IMAGES_LIST))
  IMAGING_FILE = file.path(WD, IMAGING_NAME)
  MASK_FILE = file.path(WD, MASK)
  
  if (is.null(IMAGES$V3)) {
    DATA = data.frame(SUBJECT = IMAGES$V1,
                      TP = IMAGES$V2)
    DATA = merge(DATA,
                 covars.table,
                 by = c("SUBJECT", "TP"),
                 all.x = T) %>% arrange(SUBJECT, TP)
    
  }
  else {
    DATA = data.frame(
      SUBJECT = IMAGES$V1,
      TP = IMAGES$V2,
      DEPTH = 1 - IMAGES$V3 # depth starting from cortical surface/ inverse depth
    )
    DATA = merge(DATA,
                 covars.table,
                 by = c("SUBJECT", "TP"),
                 all.x = T) %>% arrange(-DEPTH, SUBJECT, TP)
    
  }
  
  #View(DATA)
  print(DATA[!complete.cases(DATA),])
  
  imaging.mat = open_data(IMAGING_FILE, MASK_FILE)
  
  # exclude also MTX8
  #imaging.mat = imaging.mat[DATA$SYSTEM == "Classic",]
  #DATA = DATA[DATA$SYSTEM == "Classic",]
  
  excluded = which((DATA$SUBJECT %in% names(which(
    table(DATA$SUBJECT) < 4
  ))) | !complete.cases(DATA))
  
  # remove bad cases
  if (!is.null(excluded)) {
    imaging.mat = imaging.mat[-excluded,]
    DATA = DATA[-excluded,]
  }

  annotation = load_annotation(ANNOT_FILE, MASK_FILE)
  labels = unlist(read.table(LABEL_FILE))
  #labels = c('label-2','label-3','label-16','label-34','label-38','label-42','label-59','label-74', 'label-99', 'label-141')
  
  cl <- makeCluster(NPROCS)
  registerDoParallel(cl)
  accuracies = NULL
  NFOLDS = 10 #
  NREPEATS = 10
  y = unique(DATA[c("SUBJECT", "GROUP.NUM")] %>% arrange(SUBJECT))$GROUP.NUM
  
  folds = createMultiFolds(y, k = NFOLDS, times = NREPEATS)

  myperms = seq(0, NPERMS)
  iterations = expand.grid(label = labels, perm = myperms) # reorder to go for the other first
  
  k = 1
  niters = nrow(iterations)
  accuracies = matrix(NA, length(labels), length(myperms))
  colnames(accuracies) = c("value", paste0("perm_", sprintf("%03d", myperms[-1])))
  rownames(accuracies) = labels
  results = list()
  procdir = file.path(OUTDIR, 'data')
  dir.create(procdir)
  
  while (T) {
    indices = seq(NINDS * (k - 1) + 1, min(niters, NINDS * k))
    print(indices)
    if (NINDS * (k - 1) + 1 > niters)
      break
    
    get_time()
    if (T){
    train_res_list = foreach(index = indices,
                             .packages = c('caret', 'e1071', 'forecast', 'reshape', 'dplyr', 'lme4', 'brms'),
                             .export = c('prepare_data', 'open_data', 'train_model', 'get_time', 'analysis', 'Mode', 'get_pca'),
                             .combine = 'c') %dopar% analysis(iterations$label[index], imaging.mat, DATA, annotation, folds, y, HEMI,
                                                              permutate = iterations$perm[index], procdir = procdir)
    } else {
    train_res_list = NULL
    j = 1
    
    for (index in indices) {
      train_res_list[[j]] = analysis(
        iterations$label[index],
        imaging.mat,
        DATA,
        annotation,
        folds,
        y,
        HEMI,
        permutate = iterations$perm[index],
        procdir = procdir
      )
      j = j + 1
    }
    }

    #get accuracies and print out
    print("Finished run")
    accuracy_list = unlist(lapply(train_res_list, function(x)
      x$accuracy))
    for (j in seq(length(indices))) {
      index = indices[j]
      accuracies[iterations$label[index], iterations$perm[index] + 1] = accuracy_list[j]
    }
    output = format(as.data.frame(accuracies), digits = 3)
    output$label = labels
    write.table(
      output,
      file = file.path(OUTDIR, ACC_FILE),
      sep = ',',
      col.names = T,
      row.names = F,
      quote = F
    )
    
    # no need to generate every time but quick
    # binomial distribution
    binomsize = NREPEATS * length((y))
    
    output = matrix(rbinom(length(labels) * (NPERMS + 1), binomsize, p = 0.5) /
                      binomsize,
                    length(labels),
                    NPERMS + 1)
    colnames(output) = c("value", paste0("perm_", sprintf("%03d", myperms[-1])))
    output[, 'value'] = accuracies[, 'value']
    output = as.data.frame(format(output, digits = 3))
    output$label = labels
    write.table(
      output,
      file = paste0(file.path(OUTDIR, ACC_FILE), '.binom'),
      sep = ',',
      col.names = T,
      row.names = F,
      quote = F
    )
    
    results$train_res_list = c(results$train_res_list, train_res_list)
    k = k + 1
    results$accuracies = accuracies
    results$data = DATA
    OUTFILE = ifelse(is.null(HEMI), 'partial.results.rds', paste0(HEMI, '.partial.results.rds'))
    save(results, file = file.path(OUTDIR, OUTFILE))
    
  } # while
  
  stopCluster(cl)
  
  return(results)
}

