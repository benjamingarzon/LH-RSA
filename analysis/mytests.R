###############################################################
# Define the tests
###############################################################

###############################################################
# functional tests
###############################################################


testaverage = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT")
  contrast.names = c("INTERCEPT+", "INTERCEPT-")
  
  c.1        = c(1)
  c.2        = c(-1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  

  tryCatch({ 
    model = lmer(y ~ 1 + (1|SUBJECT) + (1|CONFIGURATION), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testquadraticrun = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "TRAINING.Q", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q")
  contrast.names = c("GROUP_x_TRAINING.Q+", "GROUP_x_TRAINING.Q-")
  
  c.1        = c(0, 0, 0, 0, 0, 1)
  c.2        = c(0, 0, 0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  #  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  #  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT) +  
                   (1|CONFIGURATION), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testasymptoticrun = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "TRAINING.A", "GROUP_x_TRAINING", "GROUP_x_TRAINING.A")
  contrast.names = c("GROUP_x_TRAINING.A+", "GROUP_x_TRAINING.A-")
  
  c.1        = c(0, 0, 0, 0, 0, 1)
  c.2        = c(0, 0, 0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.A) + (1 + TRAINING + TRAINING.A|SUBJECT) +  
                   (1|CONFIGURATION), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testlinearrun = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "GROUP_x_TRAINING")
  contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING-")
  
  c.1        = c(0, 0, 0, 1)
  c.2        = c(0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 

  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*(TRAINING) + (1 + TRAINING|SUBJECT) +  
                   (1|CONFIGURATION), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

modelcomparisonrun = function(y, X)
{
  
  X$y = y
  
  # if you want to select certain volumes only
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model.q = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT) +  
                   (1|CONFIGURATION), data = X)

    model.a = lmer(y ~ 1 + GROUP*(TRAINING) + (1 + TRAINING|SUBJECT) +  
                     (1|CONFIGURATION), data = X)
    
    AIC.val = AIC(model.q, model.a)
    BIC.val = BIC(model.q, model.a)    
    val = c(AIC.val$AIC[1], AIC.val$AIC[2], AIC.val$AIC[1]-AIC.val$AIC[2],
            BIC.val$BIC[1], BIC.val$BIC[2], BIC.val$BIC[1]-BIC.val$BIC[2])
    names(val) = c("AIC_quadratic", "AIC_linear", "AIC_quad-linear", "BIC_quadratic", "BIC_linear", "BIC_quad-linear")    
    return(val)
    
  },
  
  error = function(cond){
    val = rep(0, 6)
    names(val) = c("AIC_quadratic", "AIC_linear", "AIC_quad-linear", "BIC_quadratic", "BIC_linear", "BIC_quad-linear")    
    return(val)
  }
  )
  
}


###############################################################
# structural tests
###############################################################

testquadratic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "TRAINING.Q", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q")
  contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING-", "GROUP_x_TRAINING.Q+", "GROUP_x_TRAINING.Q-")
  
  c.1        = c(0, 0, 0, 0, 1, 0)
  c.2        = c(0, 0, 0, 0, -1, 0)
  c.3        = c(0, 0, 0, 0, 0, 1)
  c.4        = c(0, 0, 0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2, c.3, c.4)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  #  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  #  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
    
    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    names(val) = tags
    
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testquadratic_depth = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "DEPTH", "TRAINING", "TRAINING.Q", 
                "GROUP_x_DEPTH", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q", 
                "DEPTH_x_TRAINING", "DEPTH_x_TRAINING.Q", 
                "GROUP_x_DEPTH_x_TRAINING", "GROUP_x_DEPTH_x_TRAINING.Q")
  
  contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING.Q+", 
                     "GROUP_x_DEPTH_x_TRAINING+", "GROUP_x_DEPTH_x_TRAINING.Q+")
  
  c.1        = c(0, 0, 0, 0,   0, 0, 1, 0,   0, 0, 0, 0)
  c.2        = c(0, 0, 0, 0,   0, 0, 0, 1,   0, 0, 0, 0)
  c.3        = c(0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 1, 0)
  c.4        = c(0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 1)
  
  cont.mat = rbind(c.1, c.2, c.3, c.4)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  #  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  #  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
    
    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    names(val) = tags
    
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


modelcomparison = function(y, X)
{
  
  X$y = y
  
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model.l = lmer(y ~ 1 + GROUP*(TRAINING) + (1 + TRAINING|SUBJECT), data = X)
    
    model.a = lmer(y ~ 1 + GROUP*(TRAINING.A) + (1 + TRAINING.A|SUBJECT), data = X)
    
    model.q = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
    
    BIC.val = BIC(model.q, model.a, model.l)
    
    val = c(BIC.val$BIC[1], 
            BIC.val$BIC[2], 
            BIC.val$BIC[3])
    val = c(val, which.min(val)) # return the best model
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "Best")
    
    return(val)
  },
  
  error = function(cond){
    val = rep(0, 3)
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "Best")
    return(val)
  }
  )
  
}

modelcomparison_depth = function(y, X)
{
  
  X$y = y
  
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model.l = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING) + (1 + TRAINING|SUBJECT), data = X)
    
    model.a = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING.A) + (1 + TRAINING.A|SUBJECT), data = X)
    
    model.q = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
    
    BIC.val = BIC(model.q, model.a, model.l)
    
    val = c(BIC.val$BIC[1], 
            BIC.val$BIC[2], 
            BIC.val$BIC[3])
    val = c(val, which.min(val)) # return the best model
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "Best")
    
    return(val)
  },
  
  error = function(cond){
    val = rep(0, 3)
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "Best")
    return(val)
  }
  )
  
}




testgrouplinearhalf = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "GROUP_x_TRAINING")
  contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING-")
  
  c.1        = c(0, 0, 0, 1)
  c.2        = c(0, 0, 0, -1)
  
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
  
  # if you want to select certain volumes only
  X = subset(X, TP %in% c(1, 2, 3, 4, 5))
#  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*TRAINING + (1 + TRAINING|SUBJECT), data = X)
    
    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


test2levels = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING.L", "GROUP_x_TRAINING.L")
  contrast.names = c("GROUP_x_TRAINING.L+", "GROUP_x_TRAINING.L-")
  
  c.1        = c(0, 0, 0, 1)
  c.2        = c(0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  X = subset(X, !is.na(X$TRAINING.L))

  tryCatch({ 
#    model = lmer(y ~ 1 + GROUP*TRAINING.L + (1 + TRAINING.L|SUBJECT), data = X)
    model = lmer(y ~ 1 + GROUP*TRAINING.L + (1|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
  }
  )
}

testasymptotic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "TRAINING.A", "GROUP_x_TRAINING", "GROUP_x_TRAINING.A")
  contrast.names = c("GROUP_x_TRAINING.A+", "GROUP_x_TRAINING.A-")
  
  c.1        = c(0, 0, 0, 1, 0)
  c.2        = c(0, 0, 0, -1, 0)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + EXPERIMENT*TRAINING + TRAINING.A + (1 + TRAINING|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testanova = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("GROUP", "TRAINING", "GROUP_x_TRAINING")

  tags = c(paste0(var.names, '_F'),
           paste0(var.names, '_p'))
  X$y = y

  tryCatch({ 
    model <- aov(y ~ GROUP*TRAINING, data = X)
    Fvalues = summary(model)[[1]][ c("GROUP ", "TRAINING", "GROUP:TRAINING"), "F value"]
    pvalues = summary(model)[[1]][ c("GROUP ", "TRAINING", "GROUP:TRAINING"), "Pr(>F)"]
    val = c(Fvalues, pvalues)
    names(val) = tags
    return(val)
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    Fvalues = rep(0, length(var.names))
    val = c(Fvalues, pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

