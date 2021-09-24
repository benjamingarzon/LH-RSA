# hypothesis HF2
export_funcs = c('testlinearrun', 'testquadraticrun', 'testasymptoticrun') # export it to imagelmmr

# asymptotic tests
testgeneralization_prereg = function(y, X, ALTERNATIVE = 'greater', WAVES = seq(2, 5), select_contrasts = seq(3))
{ 
  var.names = c("INTERCEPT", "FD", "SYSTEMMTx8", "CONFIGURATION2", "CONFIGURATION3", "CONFIGURATION4",
                "GROUPExp", 
                "TRAINING", "TRAINING.A",
                "GROUPExp_x_TRAINING", "GROUPExp_x_TRAINING.A")
  
  contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_TRAINING.A-", "GROUPExp-")
  
  c.1        = c(rep(0, 9), -1, 0)
  c.2        = c(rep(0, 10), -1)
  c.3        = c(rep(0, 6), -1, rep(0, 4))

  cont.mat = rbind(c.1, c.2, c.3)
  
  cont.mat = cont.mat[select_contrasts, ]
  contrast.names = contrast.names[select_contrasts]
  
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  model = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*(TRAINING + TRAINING.A) + 
                 + (1 + TRAINING|SUBJECT), data = subset(X, WAVE %in% WAVES & CONDITION == 'UntrainedCorrect')) # + TRAINING.A
  
  tryCatch({
    
    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    
    sumglh = summary(glh, test = Ftest())
    val = c(val, sumglh$test$fstat,  sumglh$test$pvalue)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_p') 
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    tstats = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    contrast.tstats = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    val = c(val, 0, 0)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_p')
    return(val)
    
  }
  )
  
}

testfirstsessionaverage_prereg = function(y, X, ALTERNATIVE = 'greater', WAVES = seq(2, 5))
{ 
  var.names = c("INTERCEPT", "FD", "CONFIGURATION2", "CONFIGURATION3", "CONFIGURATION4")
  
  contrast.names = c("INTERCEPT")
  
  c.1        = c(1, rep(0, 4))
  
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

  model = lmer(y ~ 1 + FD + CONFIGURATION + (1 |SUBJECT), 
               data = subset(X, WAVE %in% WAVES & TP == 1)) 
  
  tryCatch({
    
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
    tstats = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    contrast.tstats = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


# asymptotic tests
testasymptoticrun = function(y, X, ALTERNATIVE = 'greater', WAVES = seq(5), select_contrasts = seq(4))
{ 
  var.names = c("INTERCEPT", "FD", "SYSTEMMTx8", "CONFIGURATION2", "CONFIGURATION3", "CONFIGURATION4",
                "GROUPExp", "CONDITIONUnt", 
                "TRAINING", "TRAINING.A",
                "GROUPExp_x_CONDITIONUnt", "GROUPExp_x_TRAINING", "GROUPExp_x_TRAINING.A",
                "CONDITIONUnt_x_TRAINING", "CONDITIONUnt_x_TRAINING.A", 
                "GROUPExp_x_CONDITIONUnt_x_TRAINING", "GROUPExp_x_CONDITIONUnt_x_TRAINING.A")
  
  contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_TRAINING.A-",
                     "GROUPExp_x_CONDITIONUnt_x_TRAINING+", "GROUPExp_x_CONDITIONUnt_x_TRAINING.A+",
                     "GROUPExp-", "GROUPExp_x_CONDITIONUnt+")
  
  c.1        = c(rep(0, 11), -1, rep(0, 5))
  c.2        = c(rep(0, 12), -1, rep(0, 4))
  c.3        = c(rep(0, 15), 1, 0)
  c.4        = c(rep(0, 16), 1)
  c.5        = c(rep(0, 6), -1, rep(0, 10))
  c.6        = c(rep(0, 10), 1, rep(0, 6))
  
  cont.mat = rbind(c.1, c.2, c.3, c.4, c.5, c.6)
  
  cont.mat = cont.mat[select_contrasts, ]
  contrast.names = contrast.names[select_contrasts]
  
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  model = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*CONDITION*(TRAINING + TRAINING.A) + 
                 + (1 + TRAINING|SUBJECT), data = subset(X, WAVE %in% WAVES)) # + TRAINING.A
  
  tryCatch({
    
    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    
    sumglh = summary(glh, test = Ftest())
    val = c(val, sumglh$test$fstat,  sumglh$test$pvalue)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_p') 
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    tstats = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    contrast.tstats = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    val = c(val, 0, 0)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_p')
    return(val)
    
  }
  )
  
}

testasymptoticrun_prereg  = function(y, X, ALTERNATIVE = 'greater') {
  return( testasymptoticrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5)))
}

testasymptoticrun_prereg_groupxtraining  = function(y, X, ALTERNATIVE = 'greater') {
  return( testasymptoticrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5), select_contrasts = c(1, 2)))
}

testasymptoticrun_prereg_groupxconditionxtraining  = function(y, X, ALTERNATIVE = 'greater') {
  return( testasymptoticrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5), select_contrasts = c(3, 4, 6)))
}

testasymptoticrun_prereg_omni  = function(y, X, ALTERNATIVE = 'greater') {
  return( testasymptoticrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5), select_contrasts = c(1, 2, 3, 4, 6)))
}



# quadratic tests
testquadraticrun = function(y, X, ALTERNATIVE = 'greater', WAVES = seq(5), select_contrasts = seq(4))
{ 
  var.names = c("INTERCEPT", "FD", "SYSTEMMTx8", "CONFIGURATION2", "CONFIGURATION3", "CONFIGURATION4",
                "GROUPExp", "CONDITIONUnt", 
                "TRAINING", "TRAINING.Q",
                "GROUPExp_x_CONDITIONUnt", "GROUPExp_x_TRAINING", "GROUPExp_x_TRAINING.Q",
                "CONDITIONUnt_x_TRAINING", "CONDITIONUnt_x_TRAINING.Q", 
                "GROUPExp_x_CONDITIONUnt_x_TRAINING", "GROUPExp_x_CONDITIONUnt_x_TRAINING.Q")
  
  contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_TRAINING.Q-",
                     "GROUPExp_x_CONDITIONUnt_x_TRAINING+", "GROUPExp_x_CONDITIONUnt_x_TRAINING.Q+",
                     "GROUPExp-", "GROUPExp_x_CONDITIONUnt+")
  
  c.1        = c(rep(0, 11), -1, rep(0, 5))
  c.2        = c(rep(0, 12), -1, rep(0, 4))
  c.3        = c(rep(0, 15), 1, 0)
  c.4        = c(rep(0, 16), 1)
  c.5        = c(rep(0, 6), -1, rep(0, 10))
  c.6        = c(rep(0, 10), 1, rep(0, 6))
  
  cont.mat = rbind(c.1, c.2, c.3, c.4, c.5, c.6)

  cont.mat = cont.mat[select_contrasts, ]
  contrast.names = contrast.names[select_contrasts]
  
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  model = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*CONDITION*(TRAINING + TRAINING.Q) + 
                 + (1 + TRAINING|SUBJECT), data = subset(X, WAVE %in% WAVES)) # + TRAINING.Q

  tryCatch({

    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    
    sumglh = summary(glh, test = Ftest())
    val = c(val, sumglh$test$fstat,  sumglh$test$pvalue)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_p') 
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    tstats = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    contrast.tstats = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
    val = c(val, 0, 0)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_p')
    return(val)
    
  }
  )
  
}

testquadraticrun_prereg  = function(y, X, ALTERNATIVE = 'greater') {
  return( testquadraticrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5)))
}

testquadraticrun_prereg_groupxtraining  = function(y, X, ALTERNATIVE = 'greater') {
  return( testquadraticrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5), select_contrasts = c(1, 2, 5)))
}

testquadraticrun_prereg_groupxconditionxtraining  = function(y, X, ALTERNATIVE = 'greater') {
  return( testquadraticrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5), select_contrasts = c(3, 4, 6)))
}

# linear tests
testlinearrun = function(y, X, ALTERNATIVE = 'greater', WAVES = seq(5), TPS = seq(7))
{ 
  
  if (max(TPS) >= 6) {
    var.names = c("INTERCEPT", "FD", "SYSTEMMTx8", 
                  "CONFIGURATION2", "CONFIGURATION3", "CONFIGURATION4",
                  "GROUPExp", "CONDITIONUnt", "TRAINING", 
                  "GROUPExp_x_CONDITIONUnt", "GROUPExp_x_TRAINING", 
                  "CONDITIONUnt_x_TRAINING", "GROUPExp_x_CONDITIONUnt_x_TRAINING")
    
    contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_CONDITIONUnt_x_TRAINING+")
    myformula = as.formula('y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*CONDITION*TRAINING + (1 + TRAINING|SUBJECT)')
    c.1        = c(rep(0, 11), -1, 0)
    c.2        = c(rep(0, 12), 1)
  } else {
    # no difference in system
    var.names = c("INTERCEPT", "FD", 
                  "CONFIGURATION2", "CONFIGURATION3", "CONFIGURATION4",
                  "GROUPExp", "CONDITIONUnt", "TRAINING", 
                  "GROUPExp_x_CONDITIONUnt", "GROUPExp_x_TRAINING", 
                  "CONDITIONUnt_x_TRAINING", "GROUPExp_x_CONDITIONUnt_x_TRAINING")
    
    contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_CONDITIONUnt_x_TRAINING+")
    myformula = as.formula('y ~ 1 + FD + CONFIGURATION + GROUP*CONDITION*TRAINING + (1 + TRAINING|SUBJECT)')
    c.1        = c(rep(0, 10), -1, 0)
    c.2        = c(rep(0, 11), 1)
  }
  
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
  model = lmer(myformula, data = subset(X, TP %in% TPS & WAVE %in% WAVES))
  
  tryCatch({
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

testlinearrun_prereg  = function(y, X, ALTERNATIVE = 'greater') {
  return( testlinearrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5)))
}

testlinearrun_prereg_half  = function(y, X, ALTERNATIVE = 'greater') {
  return( testlinearrun(y, X, ALTERNATIVE, WAVES = c(2, 3, 4, 5), TPS = seq(4)))
}

# modelcomparison
modelcomparisonrun = function(y, X)
{
  X$y = y
  X = subset(X, WAVE > 1 )
  tryCatch({ 
# not including condition CONFIGURATION
    model.l = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*TRAINING + TRAINING.Q +
                     (1 + TRAINING|SUBJECT), data = X) # + TRAINING.Q

    model.a = lmer(y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP*(TRAINING + TRAINING.A) +
                     (1 + TRAINING|SUBJECT), data = X) # + TRAINING.A
    
    model.q = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*(TRAINING + TRAINING.Q) +
                     (1 + TRAINING|SUBJECT), data = X) # + TRAINING.Q

    BIC.val = BIC(model.q, model.a, model.l)
    
    val = c(BIC.val$BIC[1], 
            BIC.val$BIC[2], 
            BIC.val$BIC[3])
    val = c(val, BIC.val$BIC[3] - BIC.val$BIC[1], which.min(val)) # return the best model
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "BIC_linear-BIC_quadratic", "Best")
    
    return(val)
    
  },
  
  error = function(cond){
    val = rep(0, 5)
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "BIC_linear-BIC_quadratic", "Best")
    return(val)
  }
  )
  
}

# testquadraticrun = function(y, X, ALTERNATIVE = 'greater')
# { 
#   var.names = c("INTERCEPT", "SYSTEM", "CONFIGURATION", "GROUP", "TRAINING", "TRAINING.Q", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q")
#   contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING.Q+")
#   
#   c.1        = c(0, 0, 0, 0, 0, 0, 1, 0)
#   c.2        = c(0, 0, 0, 0, 0, 0, 0, 1)
#   
#   cont.mat = rbind(c.1, c.2)
#   colnames(cont.mat) = var.names
#   rownames(cont.mat) = contrast.names
#   
#   tags = c(paste0(var.names, '_coef'),
#            paste0(contrast.names, '_coef'),
#            paste0(var.names, '_p'),
#            paste0(contrast.names, '_p'))
#   
#   X$y = y
# 
#   tryCatch({ 
#     model = lmer(y ~ 1  + SYSTEM + CONFIGURATION + GROUP*(TRAINING + TRAINING.Q) + 
#                     + (1 + TRAINING + TRAINING.Q|SUBJECT/TP), data = X)
#     
#     pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
#     coefs = fixef(model)
#     glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
#     contrast.pvalues = summary(glh)$test$pvalues
#     contrast.coefs = coef(glh)
#     val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
#     names(val) = tags
#     return(val)
#     
#   },
#   
#   error = function(cond){
#     # something went wrong, save special values
#     pvalues = rep(2, length(var.names))
#     coefs = rep(0, length(var.names))
#     contrast.pvalues = rep(2, length(contrast.names))
#     contrast.coefs = rep(0, length(contrast.names))
#     val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
#     names(val) = tags
#     return(val)
#     
#   }
#   )
#   
# }



testlinearrun.old = function(y, X, ALTERNATIVE = 'greater')
{ 
  
  var.names = c("INTERCEPT", "FD", "SYSTEM", "CONFIGURATION2", "CONFIGURATION3", "CONFIGURATION4", "GROUPExp", "CONDITIONUnt", "TRAINING", 
                "GROUPExp_x_CONDITIONUnt", "GROUPExp_x_TRAINING", "CONDITIONUnt_x_TRAINING", "GROUPExp_x_CONDITIONUnt_x_TRAINING")
  contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_CONDITIONUnt_x_TRAINING+")
  
  c.1        = c(rep(0, 11), -1, 0)
  c.2        = c(rep(0, 12), 1)
  
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
  model = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*CONDITION*TRAINING +
                 + (1|SUBJECT), data = X)
  
  tryCatch({
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
