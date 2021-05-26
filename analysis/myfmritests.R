# hypothesis HF2
testlinearrun = function(y, X, ALTERNATIVE = 'greater')
{ 
  
  var.names = c("INTERCEPT", "FD", "SYSTEMMTx8", "GROUPExp", "CONDITIONUnt", "TRAINING", 
                "GROUPExp_x_CONDITIONUnt", "GROUPExp_x_TRAINING", "CONDITIONUnt_x_TRAINING", "GROUPExp_x_CONDITIONUnt_x_TRAINING")
  contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_CONDITIONUnt_x_TRAINING+")
  
  c.1        = c(rep(0, 8), -1, 0)
  c.2        = c(rep(0, 9), 1)
  
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
  model = lmer(y ~ 1 + FD + SYSTEM  + GROUP*CONDITION*TRAINING +
                 + (1 + TRAINING|SUBJECT), data = X)
  
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


testquadraticrun = function(y, X, ALTERNATIVE = 'greater')
{ 
  var.names = c("INTERCEPT", "FD", "SYSTEMMTx8", "GROUPExp", "CONDITIONUnt", 
                "TRAINING", "TRAINING.Q",
                "GROUPExp_x_CONDITIONUnt", "GROUPExp_x_TRAINING", "GROUPExp_x_TRAINING.Q",
                "CONDITIONUnt_x_TRAINING", "CONDITIONUnt_x_TRAINING.Q", 
                "GROUPExp_x_CONDITIONUnt_x_TRAINING", "GROUPExp_x_CONDITIONUnt_x_TRAINING.Q")
  
  contrast.names = c("GROUPExp_x_TRAINING-", "GROUPExp_x_CONDITIONUnt_x_TRAINING+")
  
  c.1        = c(rep(0, 8), -1, rep(0, 5))
  c.2        = c(rep(0, 12), 1, 0)
  
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
  # SYSTEM
  model = lmer(y ~ 1 + FD + SYSTEM + GROUP*CONDITION*(TRAINING + TRAINING.Q) + 
                 + (1 + TRAINING + TRAINING.Q|SUBJECT), data = X)
  
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

modelcomparisonrun = function(y, X)
{
  X$y = y
  
  tryCatch({ 
# not including condition CONFIGURATION
    model.l = lmer(y ~ 1 + FD + SYSTEM  + GROUP*TRAINING + TRAINING.Q +
                     (1 + TRAINING|SUBJECT), data = X)

    model.a = lmer(y ~ 1 + FD + SYSTEM  + GROUP*(TRAINING + TRAINING.A) +
                     (1 + TRAINING |SUBJECT), data = X)
    
    model.q = lmer(y ~ 1 + FD + SYSTEM + GROUP*(TRAINING + TRAINING.Q) +
                     (1 + TRAINING|SUBJECT), data = X)

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
