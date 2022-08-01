# hypothesis HF2
export_funcs = c('testquadraticrun',
                 'testasymptoticrun',
                 'testcubicrun',
                 'modelcomparisonrun') # export it to imagelmmr


testasymptoticrun = function(y, X, WAVES = seq(5))
{
  var.names = c(
    "FD",
    "SYSTEMMTx8",
    "CONFIGURATION",
    "GROUP",
    "CONDITION",
    "TRAINING.A",
    "TRAINING",
    "GROUP_x_CONDITION",
    "GROUP_x_TRAINING.A",
    "CONDITION_x_TRAINING.A",
    "GROUP_x_CONDITION_x_TRAINING.A"
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
    y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP * CONDITION * TRAINING.A + poly(TRAINING, 3, raw = T) +
      (1 + TRAINING | SUBJECT),
    data = X, contrasts=list(GROUP=contr.sum, CONDITION=contr.sum)
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


testquadraticrun = function(y, X, WAVES = seq(5))
{
  var.names = c(
    "FD",
    "SYSTEMMTx8",
    "CONFIGURATION",
    "GROUP",
    "CONDITION",
    "TRAINING.Q",
    "TRAINING",
    "GROUP_x_CONDITION",
    "GROUP_x_TRAINING.Q",
    "CONDITION_x_TRAINING.Q",
    "GROUP_x_CONDITION_x_TRAINING.Q"
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
    y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP * CONDITION * TRAINING.Q + poly(TRAINING, 3, raw = T) +
      (1 + TRAINING | SUBJECT),
    data = X, contrasts=list(GROUP=contr.sum, CONDITION=contr.sum)
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


testcubicrun = function(y, X, WAVES = seq(5))
{
  var.names = c(
    "FD",
    "SYSTEMMTx8",
    "CONFIGURATION",
    "GROUP",
    "CONDITION",
    "TRAINING",
    "GROUP_x_CONDITION",
    "GROUP_x_TRAINING",
    "CONDITION_x_TRAINING",
    "GROUP_x_CONDITION_x_TRAINING"
  )
  
  tags = c(paste0(var.names, '_chisq'),
           paste0(var.names, '_p'))

  X$y = y
  
  X = subset(X, WAVE %in% WAVES)
  
  X$TRAINING = scale((X$TP - 1)*6, center = F, scale = T)
  
  model = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * CONDITION * poly(TRAINING, 3, raw = T) + 
      (1 + TRAINING | SUBJECT), data = X, contrasts=list(GROUP=contr.sum, CONDITION=contr.sum))
  
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


modelcomparisonrun = function(y, X, WAVES = seq(5))
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
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * CONDITION * TRAINING + poly(TRAINING, deg, raw = T) +
        (1 + TRAINING | SUBJECT),
      data = X
    )
    
    model.a = lmer(
      y ~ 1 + FD + SYSTEM  + CONFIGURATION + GROUP * CONDITION * TRAINING.A + poly(TRAINING, deg, raw = T) +
        (1 + TRAINING | SUBJECT),
      data = X
    )
    
    model.q = lmer(
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * CONDITION * TRAINING.Q + poly(TRAINING, deg, raw = T) +
        (1 + TRAINING | SUBJECT),
      data = X
    )
    
    model.q2 = lmer(
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * CONDITION * poly(TRAINING, 2, raw = T) +  poly(TRAINING, deg, raw = T) +
        (1 + TRAINING | SUBJECT),
      data = X
    )
    
    model.q3 = lmer(
      y ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP * CONDITION * poly(TRAINING, deg, raw = T) +
        (1 + TRAINING | SUBJECT),
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


testasymptoticrun_prereg  = function(y, X) {
  return(testasymptoticrun(y, X, WAVES = c(2, 3, 4, 5)))
}


testquadraticrun_prereg  = function(y, X) {
  return(testquadraticrun(y, X, WAVES = c(2, 3, 4, 5)))
}

testcubicrun_prereg  = function(y, X) {
  return(testcubicrun(y, X, WAVES = c(2, 3, 4, 5)))
}


modelcomparisonrun_prereg  = function(y, X) {
  return(modelcomparisonrun(y, X, WAVES = c(2, 3, 4, 5)))
}

