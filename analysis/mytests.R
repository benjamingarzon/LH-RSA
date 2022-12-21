###############################################################
# Define the tests
###############################################################

###############################################################
# structural tests
###############################################################

# structural hypotheses

testasymptotic = function(y, X, WAVES = seq(5))
{
  var.names = c(
    "INTERCEPT",
    "FD",
    "SYSTEMMTx8",
    "CONFIGURATION",
    "GROUP",
    "TRAINING.A",
    "TRAINING.C",
    "GROUP_x_TRAINING.A"
  )
  
  
  tags = c(paste0(var.names, '_chisq'),
           paste0(var.names, '_p'))
  
  training = seq(0, 6) * 6
  training.quadratic = -(training - max(training) / 2) ^ 2 # center around midpoint
  min.quad = min(training.quadratic)
  training.quadratic = training.quadratic - min.quad
  training.asymptotic = c(training.quadratic[1:4], rep(training.quadratic[4], 3))
  
  X$y = y
  
  X = subset(X, WAVE %in% WAVES)
  
  X$TRAINING = scale(training[X$TP], center = F, scale = T)
  X$TRAINING.A = scale(training.asymptotic[X$TP], center = F, scale = T)
  
  model = lmer(
    y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP * TRAINING.A + poly(TRAINING, 3, raw = T) +
      (1 | SUBJECT),
    data = X, contrasts=list(GROUP=contr.sum)
  )
  
  tryCatch({
    an = Anova(model, type = "III")
    pvalues = unlist(an[, "Pr(>Chisq)"])
    stats = unlist(an[, "Chisq"])
    val = c(stats, pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond) {
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    stats = rep(0, length(var.names))
    val = c(stats, pvalues)
    names(val) = tags
    return(val)
    
  })
}


testquadratic = function(y, X, WAVES = seq(5))
{
  var.names = c(
    "INTERCEPT",
    "FD",
    "SYSTEMMTx8",
    "CONFIGURATION",
    "GROUP",
    "TRAINING.Q",
    "TRAINING.C",
    "GROUP_x_TRAINING.Q"
  )
  
  
  tags = c(paste0(var.names, '_chisq'),
           paste0(var.names, '_p'))
  
  training = seq(0, 6) * 6
  training.quadratic = -(training - max(training) / 2) ^ 2 # center around midpoint
  min.quad = min(training.quadratic)
  training.quadratic = training.quadratic - min.quad
  
  X$y = y
  
  X = subset(X, WAVE %in% WAVES)
  
  X$TRAINING = scale(training[X$TP], center = F, scale = T)
  X$TRAINING.Q = scale(training.quadratic[X$TP], center = F, scale = T)
  
  model = lmer(
    y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP * TRAINING.Q + poly(TRAINING, 3, raw = T) +
      (1 | SUBJECT),
    data = X, contrasts=list(GROUP=contr.sum)
  )
  
  tryCatch({
    an = Anova(model, type = "III")
    pvalues = unlist(an[, "Pr(>Chisq)"])
    stats = unlist(an[, "Chisq"])
    val = c(stats, pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond) {
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    stats = rep(0, length(var.names))
    val = c(stats, pvalues)
    names(val) = tags
    return(val)
    
  })
}

testlinear = function(y, X, WAVES = seq(5))
{
  var.names = c(
    "INTERCEPT",
    "FD",
    "SYSTEMMTx8",
    "CONFIGURATION",
    "GROUP",
    "TRAINING",
    "TRAINING.C",
    "GROUP_x_TRAINING"
  )
  
  
  tags = c(paste0(var.names, '_chisq'),
           paste0(var.names, '_p'))
  
  training = seq(0, 6) * 6

  X$y = y
  
  X = subset(X, WAVE %in% WAVES)
  
  X$TRAINING = scale(training[X$TP], center = F, scale = T)

  model = lmer(
    y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP * TRAINING + poly(TRAINING, 3, raw = T) +
      (1 | SUBJECT),
    data = X, contrasts=list(GROUP=contr.sum)
  )
  
  tryCatch({
    an = Anova(model, type = "III")
    pvalues = unlist(an[, "Pr(>Chisq)"])
    stats = unlist(an[, "Chisq"])
    val = c(stats, pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond) {
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    stats = rep(0, length(var.names))
    val = c(stats, pvalues)
    names(val) = tags
    return(val)
    
  })
}

modelcomparison = function(y, X, WAVES = seq(5))
{
  training = seq(0, 6) * 6
  training.quadratic = -(training - max(training) / 2) ^ 2 # center around midpoint
  min.quad = min(training.quadratic)
  training.quadratic = training.quadratic - min.quad
  training.asymptotic = c(training.quadratic[1:4], rep(training.quadratic[4], 3))
  
  X$y = y
  
  X = subset(X, WAVE %in% WAVES)
  
  X$TRAINING = scale(training[X$TP], center = F, scale = T)
  X$TRAINING.Q = scale(training.quadratic[X$TP], center = F, scale = T)
  X$TRAINING.A = scale(training.asymptotic[X$TP], center = F, scale = T)
  
  deg = 3
  tryCatch({
    model.l = lmer(
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * TRAINING + poly(TRAINING, deg, raw = T) +
        (1 | SUBJECT),
      data = X
    )
    
    model.a = lmer(
      y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP * TRAINING.A + poly(TRAINING, deg, raw = T) +
        (1 | SUBJECT),
      data = X
    )
    
    model.q = lmer(
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * TRAINING.Q + poly(TRAINING, deg, raw = T) +
        (1 | SUBJECT),
      data = X
    )
    
    model.q2 = lmer(
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * poly(TRAINING, 2, raw = T) +  poly(TRAINING, deg, raw = T) +
        (1 | SUBJECT),
      data = X
    )
    
    model.q3 = lmer(
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * poly(TRAINING, deg, raw = T) +
        (1 | SUBJECT),
      data = X
    )
    
    BIC.val = BIC(model.l, model.a, model.q, model.q2, model.q3)
    
    val = c(BIC.val$BIC[1],
            BIC.val$BIC[2],
            BIC.val$BIC[3],
            BIC.val$BIC[4],
            BIC.val$BIC[5])
    
    val = c(val, BIC.val$BIC[2] - BIC.val$BIC[3], which.min(val)) # return the best model
    names(val) = c(
      "BIC_linear",
      "BIC_asymptotic",
      "BIC_quadratic",
      "BIC_quadratic2",
      "BIC_quadratic3",
      "BIC_asymptotic-BIC_quadratic",
      "Best"
    )
    
    return(val)
    
  },
  
  error = function(cond) {
    val = rep(0, 7)
    names(val) = c(
      "BIC_linear",
      "BIC_asymptotic",
      "BIC_quadratic",
      "BIC_quadratic2",
      "BIC_quadratic3",
      "BIC_asymptotic-BIC_quadratic",
      "Best"
    )
    return(val)
  })
  
}

reliability = function(y, X)
{
  
  X$y = y
  X = X %>% filter(GROUP == 'Control' & SYSTEM == 'Classic') 
  tryCatch({ 
    
    XX = cast(X, 'SUBJECT ~ TP', value = 'y')
    myICC = icc(XX[, -1], model = "twoway", type = "agreement")
    val = myICC$value
    names(val) = c("ICC")
    
    return(val)
    
  },
  
  error = function(cond){
    val = 0
    names(val) = c("ICC")
    return(val)
  }
  )
  
}

reliability_depth = function(y, X)
{
  
  X$y = y
  X = X %>% filter(GROUP == 'Control' & SYSTEM == 'Classic' & DEPTH < 1) 
  tryCatch({ 
    depths = unique(X$DEPTH)
    vals = rep(0, length(depths))
    names(vals) = depths
    for (depth in depths){
      print(depth)
     XX = cast(X %>% filter(DEPTH == depth), 'SUBJECT ~ TP', value = 'y')
     myICC = icc(XX[, -1], model = "twoway", type = "agreement")
     vals[as.character(depth)] = myICC$value
    }
    val = mean(vals, na.rm =T)
    names(val) = c("ICC")
    
    return(val)
    
  },
  
  error = function(cond){
    val = 0
    names(val) = c("ICC")
    return(val)
  }
  )
  
}


## same tests, including depth to study T1
modelcomparison_depth = function(y, X)
{
  
  X$y = y
  X = subset(X, DEPTH < 1)
  tryCatch({ 
    model.l3 = lmer(y ~ 1 + FD + SYSTEM + DEPTH + GROUP*TRAINING + TRAINING.A  + (1 |SUBJECT), data = X)
    model.l2 = lmer(y ~ 1 + FD + SYSTEM + DEPTH + GROUP*TRAINING + (1 |SUBJECT), data = X)
    model.l = lmer(y ~ 1 + FD + SYSTEM + DEPTH + GROUP*TRAINING + TRAINING.Q  + (1 |SUBJECT), data = X)
    model.a = lmer(y ~ 1 + FD + SYSTEM + DEPTH + GROUP*(TRAINING + TRAINING.A) + (1 |SUBJECT), data = X)
    model.q = lmer(y ~ 1 + FD + SYSTEM + DEPTH + GROUP*(TRAINING + TRAINING.Q) + (1 |SUBJECT), data = X)
    
    BIC.val = BIC(model.q, model.a, model.l, model.l2, model.l3)
    
    val = c(BIC.val$BIC[1], 
            BIC.val$BIC[2], 
            BIC.val$BIC[3],
            BIC.val$BIC[4],
            BIC.val$BIC[5])
  val = c(val, BIC.val$BIC[3] - BIC.val$BIC[1], which.min(val)) # return the best model
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "BIC_linear2", "BIC_linear3", "BIC_linear-BIC_quadratic", "Best")
    
    return(val)
  },
  
  error = function(cond){
    val = rep(0, 7)
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "BIC_linear2", "BIC_linear3", "BIC_linear-BIC_quadratic", "Best")
    return(val)
  }
  )
  
}

testquadratic_depth = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "FD", "SYSTEM", "DEPTH", "GROUP" , "TRAINING", "TRAINING.Q", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q")
  
  contrast.names = c(
    "GROUP_x_TRAINING+", 
    "GROUP_x_TRAINING.Q+" 
  )
  
  c.1        = c(0, 0, 0, 0, 0, 0, 0,  1,  0)
  c.2        = c(0, 0, 0, 0, 0, 0,  0,  0,  1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  X = subset(X, DEPTH < 1)
  
  tryCatch({ 
    model = lmer(y ~ 1 + FD + SYSTEM + DEPTH + GROUP*(TRAINING + TRAINING.Q) + (1 |SUBJECT), data = X)
    model_type = 1

    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE)
    contrast.pvalues = as.numeric(summary(glh)$test$pvalues)
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    names(contrast.pvalues) = names(contrast.coefs)
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    
    sumglh = summary(glh, test = Ftest())
    val = c(val, sumglh$test$fstat, sumglh$test$SSH, sumglh$test$pvalue, model_type)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_coef', 'Omni_p', 'modeltype_coef')
    
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    val = c(val, rep(0, 3), -1)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_coef', 'Omni_p', 'modeltype_coef')
    return(val)
    
  }
  )
  
}

testlinear_depth = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "FD", "SYSTEM", "DEPTH", "GROUP" , "TRAINING", "GROUP_x_TRAINING")
  
  contrast.names = c("GROUP_x_TRAINING+")
  
  c.1        = c(rep(0, 6), 1)
  
  cont.mat = rbind(c.1)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  X = subset(X, DEPTH < 1)
  
  tryCatch({ 
    
    model = lmer(y ~ 1 + FD + SYSTEM + DEPTH + GROUP*TRAINING + (1 |SUBJECT), data = X)
    model_type = 1

    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE)
    contrast.pvalues = as.numeric(summary(glh)$test$pvalues)
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    names(contrast.pvalues) = names(contrast.coefs)
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    
    sumglh = summary(glh, test = Ftest())
    val = c(val, sumglh$test$fstat, sumglh$test$SSH, sumglh$test$pvalue, model_type)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_coef', 'Omni_p', 'modeltype_coef')
    
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    val = c(val, rep(0, 3), -1)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_coef', 'Omni_p', 'modeltype_coef')
    return(val)
    
  }
  )
  
}

modelcomparison_prereg  = function(y, X) {
  return(modelcomparison(y, X, WAVES = c(2, 3, 4, 5)))
}

testasymptotic_prereg  = function(y, X) {
  return(testasymptotic(y, X, WAVES = c(2, 3, 4, 5)))
}

testquadratic_prereg  = function(y, X) {
  return(testquadratic(y, X, WAVES = c(2, 3, 4, 5)))
}

testlinear_prereg  = function(y, X) {
  return(testlinear(y, X, WAVES = c(2, 3, 4, 5)))
}

