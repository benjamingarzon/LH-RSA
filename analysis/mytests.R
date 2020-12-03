###############################################################
# Define the tests
###############################################################
###############################################################
# functional tests
###############################################################
testaverage = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "CONFIGURATION")
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
    model = lmer(y ~ 1 + CONFIGURATION + (1|SUBJECT), data = X)
    
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

# hypothesis HF2
testlinearrun = function(y, X, ALTERNATIVE = 'greater')
{ 
   var.names = c("INTERCEPT", "SYSTEM", "CONFIGURATION", "GROUP", "TRAINING", "GROUP_x_TRAINING")
   contrast.names = c("GROUP_x_TRAINING-")

  c.1        = c(0, 0, 0,  0, 0, -1)

  cont.mat = rbind(c.1)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names

  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))

  X$y = y

  tryCatch({
    model = lmer(y ~ 1 + SYSTEM + CONFIGURATION + GROUP*TRAINING +
                    + (1 + TRAINING |SUBJECT/TP), data = X)

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
  
  tryCatch({ 
    model.q = lmer(y ~ 1 + SYSTEM + CONFIGURATION + GROUP*(TRAINING + TRAINING.Q) + 
                    + (1 + TRAINING + TRAINING.Q|SUBJECT/TP), data = X)

    model.a = lmer(y ~ 1 + SYSTEM + CONFIGURATION + GROUP*(TRAINING + TRAINING.A) + 
                     + (1 + TRAINING + TRAINING.A|SUBJECT/TP), data = X)
    
    model.l = lmer(y ~ 1 + SYSTEM + CONFIGURATION + GROUP*TRAINING + 
                      + (1 + TRAINING|SUBJECT/TP), data = X)
    
    BIC.val = BIC(model.q, model.a, model.l)
    
    val = c(BIC.val$BIC[1], 
            BIC.val$BIC[2], 
            BIC.val$BIC[3])
    val = c(val, which.min(val)) # return the best model
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "Best")
    
    
    return(val)
    
  },
  
  error = function(cond){
    val = rep(0, 4)
    names(val) = c("BIC_quadratic", "BIC_asymptotic", "BIC_linear", "Best")
    return(val)
  }
  )
  
}


###############################################################
# multivariate indices
###############################################################

# hypothesis HF1
testlinear_mv = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "SYSTEM", "CONFIGURATION", "GROUP", "TRAINING", "GROUP_x_TRAINING")
  contrast.names = c("GROUP_x_TRAINING-")
  
  c.1        = c(0, 0, 0,  0, 0, -1)

  cont.mat = rbind(c.1)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  X$y = y
  
  tryCatch({ 
    model = lmer(y ~ 1 + SYSTEM + CONFIGURATION + GROUP*TRAINING +
                   + (1 + TRAINING |SUBJECT), data = X)
    

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

###############################################################
# structural tests
###############################################################
testvariance = function(y, X, ALTERNATIVE = 'greater')
{
  
  X$y = y
  
  X.subject = X %>% filter(SYSTEM == 'Classic') %>% group_by(SUBJECT, GROUP) %>% dplyr::summarise(y.sd = sd(y))

  tags = c('sddiff_coef', 'sddiff_tstat', 'sddiff_p')
  
  tryCatch({ 
    model.t = t.test(y.sd ~ as.factor(GROUP), data = X.subject, alternative = ALTERNATIVE)
    
    pvalue = model.t$p.value
    coef = diff(model.t$estimate)
    tstat = model.t$statistic
    
    val = c(coef, tstat, pvalue)
    names(val) = tags
    return(val)
  },
  
  error = function(cond){
    val = rep(0, 3)
    names(val) = tags
    return(val)
  }
  )
  
}


# hypotheses HS1.1 and HS1.2
testlinear = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "SYSTEM", "GROUP" , "TRAINING", "GROUP_x_TRAINING")

  contrast.names = c("GROUP_x_TRAINING+")
  
  # contrast.names = c(
  #   "GROUP+", 
  #   "GROUP-",
  #   "GROUP_x_TRAINING+", 
  #   "GROUP_x_TRAINING-")
  
#  c.1        = c(0, 0, 1, 0,  0)
#  c.2        = c(0, 0,-1, 0,  0)
  c.3        = c(0, 0, 0, 0,  1)
#  c.4        = c(0, 0, 0, 0, -1)
  
#  cont.mat = rbind(c.1, c.2, c.3, c.4)
  cont.mat = rbind(c.3)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  tryCatch({ 
#    model = lmer(y ~ 1 + SYSTEM + GROUP*TRAINING + (1 + TRAINING |SUBJECT), data = X)
    model = lmer(y ~ 1 + SYSTEM + GROUP*TRAINING + (1 |SUBJECT), data = X)
    
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

# hypotheses HS1.1 and HS1.2
testlinear_Classic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP" , "TRAINING", "GROUP_x_TRAINING")
  
  contrast.names = c("GROUP_x_TRAINING+")
  
  c.1       = c(0, 0, 0, 0,  1)

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
  X = X %>% filter(SYSTEM == "Classic") %>% group_by(SUBJECT) %>% mutate(n = sum(!is.na(y))) %>% filter(n>4)
  tryCatch({ 
    model = lmer(y ~ 1 +  GROUP*TRAINING + (1 |SUBJECT), data = X)
    
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

testlinear_Classic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP" , "TRAINING", "GROUP_x_TRAINING")
  
  contrast.names = c("GROUP_x_TRAINING+")
  
  c.1       = c(0, 0, 0, 0,  1)
  
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
  X = X %>% group_by(SUBJECT) %>% mutate(n = sum(!is.na(y)))
  tryCatch({ 
    X = subset(X, SYSTEM == "Classic" & n > 4)
    model = lmer(y ~ 1 +  GROUP*TRAINING + (1 |SUBJECT), data = X)
    
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


# hypotheses HS1.1 and HS1.2

testquadratic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "SYSTEM", "GROUP" , "TRAINING", "TRAINING.Q", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q")

  contrast.names = c(
    #"GROUP_x_TRAINING+", 
    "GROUP_x_TRAINING-", 
    #"GROUP_x_TRAINING.Q+" 
    "GROUP_x_TRAINING.Q-"
    )

  #c.1        = c(0, 0, 0, 0, 0,  1,  0)
  c.2        = c(0, 0, 0, 0, 0,  -1,  0)
#  c.3        = c(0, 0, 0, 0,  0,  0,  1)
  c.4        = c(0, 0, 0, 0,  0, 0,  -1)
  
  cont.mat = rbind(#c.1, 
                   c.2, 
                   #c.3#, 
                   c.4
                   )
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  tryCatch({ 
#    model = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
    model = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.Q) + (1|SUBJECT), data = X)
    
    coefs = fixef(model)
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    tstats = summary(model)$coefficients[ , "t value"]
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE)
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    contrast.tstats = summary(glh)$test$tstat
    val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
      
    sumglh = summary(glh, test = Ftest())
    val = c(val, sumglh$test$fstat, sumglh$test$coef, sumglh$test$pvalue)
    names(val) = c(tags, 'OmniF_tstat', 'Omni_coef', 'Omni_p')
    
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    val = c(val, 0, 0)
    names(val) = c(tags, 'Omni_F', 'Omni_p')
    return(val)
    
  }
  )
  
}

testquadratic_Classic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP" , "TRAINING", "TRAINING.Q", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q")
  
  contrast.names = c(
    "GROUP_x_TRAINING-", 
    "GROUP_x_TRAINING.Q-"
  )
  
  c.1        = c(0, 0, 0, 0, -1,  0)
  c.3        = c(0, 0, 0, 0,  0, -1)
  
  cont.mat = rbind(c.1, 
                   c.3
  )
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  X = X %>% filter(SYSTEM == "Classic") %>% group_by(SUBJECT) %>% mutate(n = sum(!is.na(y))) %>% filter(n>4)
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1|SUBJECT), data = X)
    
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

# hypotheses HS1.3 and HS2
testquadratic_depth = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "SYSTEM", "GROUP", "DEPTH", "TRAINING", "TRAINING.Q", 
                "GROUP_x_DEPTH", "GROUP_x_TRAINING", "GROUP_x_TRAINING.Q", 
                "DEPTH_x_TRAINING", "DEPTH_x_TRAINING.Q", 
                "GROUP_x_DEPTH_x_TRAINING", "GROUP_x_DEPTH_x_TRAINING.Q")
  
  contrast.names = c("GROUP_x_TRAINING-", "GROUP_x_TRAINING.Q-")#, 
#                     "GROUP_x_DEPTH_x_TRAINING-", "GROUP_x_DEPTH_x_TRAINING.Q-")
  
  c.1        = c(0, 0, 0, 0, 0, 0,   0, -1, 0,   0, 0,  0, 0)
  c.2        = c(0, 0, 0, 0, 0, 0,   0, 0, -1,   0, 0,  0, 0)
#  c.3        = c(0, 0, 0, 0, 0, 0,   0, 0,  0,   0, 0,  -1, 0)
#  c.4        = c(0, 0, 0, 0, 0, 0,   0, 0,  0,   0, 0,  0, -1)
  
  cont.mat = rbind(c.1, c.2) #, c.3, c.4)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_tstat'),
           paste0(contrast.names, '_tstat'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  tryCatch({ 
  #  model = lmer(y ~ 1 + SYSTEM + GROUP*DEPTH*(TRAINING + TRAINING.Q)  + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
    model = lmer(y ~ 1 + SYSTEM + GROUP*DEPTH*(TRAINING + TRAINING.Q)  + (1|SUBJECT), data = X)
    
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


modelcomparison = function(y, X)
{
  
  X$y = y
  
  #X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  #X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  #X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
#    model.l = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING)  + (1 + TRAINING|SUBJECT), data = X)
#    model.a = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.A) + (1 + TRAINING + TRAINING.A|SUBJECT), data = X)
#    model.q = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
    model.l = lmer(y ~ 1 + SYSTEM + GROUP*TRAINING + TRAINING.Q  + (1 |SUBJECT), data = X)
    model.a = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.A) + (1 |SUBJECT), data = X)
    model.q = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.Q) + (1 |SUBJECT), data = X)
    
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

modelcomparison_depth = function(y, X)
{
  
  X$y = y
  
  #X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  #X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  #X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    # model.l = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING) + SYSTEM + (1 + TRAINING|SUBJECT), data = X)
    # model.a = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.A) + SYSTEM  + (1 + TRAINING + TRAINING.A|SUBJECT), data = X)
    # model.q = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.Q) + SYSTEM  + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)

    model.l = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING) + DEPTH*TRAINING.Q + SYSTEM + (1|SUBJECT), data = X)
    model.a = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.A) + SYSTEM  + (1|SUBJECT), data = X)
    model.q = lmer(y ~ 1 + GROUP*DEPTH*(TRAINING + TRAINING.Q) + SYSTEM  + (1|SUBJECT), data = X)
    
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


###############################################################
# Other
###############################################################

# Test for an effect of receiver system
testsystem = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "SYSTEM")
  contrast.names = c("SYSTEM+", "SYSTEM-")
  c.1        = c(0, 1)
  c.2        = c(0, -1)

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
  
  # test only in subjects with more than one system
  X = subset(X, SUBJECT %in% MTX8_subjects)  
  
  tryCatch({ 
    model = lmer(y ~ 1 + SYSTEM + (1|SUBJECT), data = X)
    
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


# Associations with IQ

testreasoning = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "REASONING", "GENDER", "SYSTEM", "LIQUID") 
  contrast.names = c("REASONING+", "REASONING-")
  
  c.1        = c(0,  1, 0, 0, 0)
  c.2        = c(0, -1, 0, 0, 0)

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
  
  tryCatch({ 
    model = lmer(y ~ 1 + REASONING + GENDER + SYSTEM + LIQUID + (1|SUBJECT), data = X) # SYSTEM
    
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

testreasoning_depth = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "REASONING", "DEPTH", "REASONING_x_DEPTH", "LIQUID", "SYSTEM")
  contrast.names = c("REASONING+", "REASONING-", "REASONING_x_DEPTH+", "REASONING_x_DEPTH-")
  
  c.1        = c(0,  1, 0, 0, 0)
  c.2        = c(0, -1, 0, 0, 0)
  
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
  
  
  tryCatch({ 
    model = lmer(y ~ 1 + REASONING*DEPTH + LIQUID + SYSTEM + (1|SUBJECT), data = X)
    
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


testgrouplinearhalf = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "SYSTEM", "GROUP", "TRAINING", "GROUP_x_TRAINING")
  contrast.names = c("GROUP_x_TRAINING+")
  
  c.1        = c(0, 0, 0, 0, 1)

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
  
  # if you want to select certain volumes only
  X = subset(X, TP %in% c(1, 2, 3, 4, 5))
  #  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + SYSTEM + GROUP*TRAINING + (1 |SUBJECT), data = X)
    
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


testgrouplinearhalf_depth = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "DEPTH", "SYSTEM", "GROUP", "TRAINING", "GROUP_x_TRAINING")
  contrast.names = c("GROUP_x_TRAINING-")
  
  c.1        = c(0, 0, 0, 0, 0, -1)
  
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
  
  # if you want to select certain volumes only
  X = subset(X, TP %in% c(1, 2, 3, 4, 5))
  #  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + DEPTH + SYSTEM + GROUP*TRAINING + (1 + TRAINING|SUBJECT), data = X)
    
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

###############################################################
# not used
###############################################################


# testgrouplinearhalf = function(y, X, ALTERNATIVE = 'greater')
# {
#   var.names = c("INTERCEPT", "GROUP", "TRAINING", "GROUP_x_TRAINING")
#   contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING-")
#   
#   c.1        = c(0, 0, 0, 1)
#   c.2        = c(0, 0, 0, -1)
#   
#   cont.mat = rbind(c.1, c.2)
#   colnames(cont.mat) = var.names
#   rownames(cont.mat) = contrast.names
#   
#   tags = c(paste0(var.names, '_coef'),
#            paste0(contrast.names, '_coef'),
#            paste0(var.names, '_tstat'),
#            paste0(contrast.names, '_tstat'),
#            paste0(var.names, '_p'),
#            paste0(contrast.names, '_p'))
#   
#   X$y = y
#   
#   # if you want to select certain volumes only
#   X = subset(X, TP %in% c(1, 2, 3, 4, 5))
# #  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
#   
#   tryCatch({ 
#     model = lmer(y ~ 1 + GROUP*TRAINING + (1 + TRAINING|SUBJECT), data = X)
#     
#     coefs = fixef(model)
#     pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
#     tstats = summary(model)$coefficients[ , "t value"]
#     glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
#     contrast.pvalues = summary(glh)$test$pvalues
#     contrast.coefs = coef(glh)
#     contrast.tstats = summary(glh)$test$tstat
#     val = c(coefs, contrast.coefs, tstats, contrast.tstats, pvalues, contrast.pvalues)
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

