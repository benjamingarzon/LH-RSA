library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(lme4)
library(lmerTest)
library(mgcv)
STARSIZE = 5
k.SESS = 7
k.CONFIG = 7
smooth = "tp"
clean_same = F
contrast = NULL
# correct for FD. FD vs type of 
# correct for number of correct trials

do_tests = function(mydata, par, myfile, odd, analysis_type = 'groupxtrainingxmeasure', output_text = F, use_GAM = T){
  
  unlink(myfile)
  mydata = mydata %>% filter(! MEASURE %in% odd)
  labels = sort(unique(mydata$label))
  groups = unique(mydata$GROUP)
  #  browser()
  #mydata$MEASURE = relevel(factor(mydata$MEASURE), ref = 'Untrained')
  
  infos = NULL
  mylines = c(par)
  if (!use_GAM){
    for (l in labels) {
      if (analysis_type == 'groupxtrainingxmeasure') 

        model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + (GROUP*TRAINING*MEASURE) + (1 |subject), data = mydata %>% 
                       filter(label == l))
        if (analysis_type == 'groupxtraining') 
        model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + (GROUP*TRAINING) + (1 |subject), data = mydata %>% 
                       filter(label == l))
      if (analysis_type == 'training') 
        model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + TRAINING + (1|subject), data = mydata %>% 
                       filter(label == l))
      if (analysis_type == 'groupxmeasure') 
        model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*MEASURE + (1 |subject), data = mydata %>% 
                       filter(label == l))
      if(isSingular(model)) {
        browser()
        print(myfile)
        print(summary(model))
        
      }
      
      infos = rbind(infos, summary(model)$coefficients[par, ])
    }
    
  } else {
    # use GAM
    for (l in labels) {
      #print(l)
      if (analysis_type == 'groupxtrainingxmeasure') 
        model = gam(value ~ 1 + FD + SYSTEM + CONFIGURATION + (GROUP*TRAINING*MEASURE) +
                      s(subject, bs = "re") +
                      s(TRAINING, by = interaction(GROUP, MEASURE), k = k.SESS, bs = smooth),
                    data = mydata %>% filter(label == l) %>% mutate(subject = factor(subject)), 
                    method = 'REML')
      
      if (analysis_type == 'groupxmeasure') 
        model = gam(value ~ 1 + FD + SYSTEM + CONFIGURATION + (GROUP*MEASURE) +
                      s(subject, bs = "re") +
                      s(TRAINING, by = interaction(GROUP, MEASURE), k = k.SESS, bs = smooth),
                    data = mydata %>% filter(label == l) %>% mutate(subject = factor(subject)), 
                    method = 'REML')
      #gam.check(model)
      
      if (analysis_type == 'groupxtraining') 
        model = gam(value ~ 1 + FD + SYSTEM + CONFIGURATION + (GROUP*TRAINING) +
                      s(subject, bs = "re") +
                      s(TRAINING, by = interaction(GROUP), k = k.SESS, bs = smooth), 
                    data = mydata %>% filter(label == l) %>% mutate(subject = factor(subject)), 
                    method = 'REML')
      
      if (analysis_type == 'training') 
        model = gam(value ~ 1 + FD + SYSTEM + CONFIGURATION + TRAINING +
                      s(subject, bs = "re") +
                      s(TRAINING, k = k.SESS, bs = smooth), 
                    data = mydata %>% filter(label == l) %>% mutate(subject = factor(subject)), 
                    method = 'REML')
      #print(summary(model))
      mysum = summary(model, freq = T)
      
      p.table = mysum$p.table
      p.table = cbind(p.table, df = df.residual(model))
      p.table = p.table[, c('Estimate', 'Std. Error', 'df', 't value', 'Pr(>|t|)')]
      print(par)
      print(rownames(p.table))
      infos = rbind(infos, p.table[par, ])
    }
    
    #print(
    #  ggplot(unpaced_trials.wrong, aes(x = sess_num, y = wrong_trials.pred, col = group)) + geom_point() + geom_line() +
    #    facet_grid(CONFIGURATION.SIMPLE ~ seq_train)
    #)
    
  }
  
  infos = cbind(infos, p.adjust(infos[, 'Pr(>|t|)'], method = 'fdr'))
  significance = ifelse(infos[, 5] < 0.05, '*', 'n. s.')
  significance[infos[, 6] < 0.05] = '**'
  infos = cbind(infos, significance)
  colnames(infos)[6:7] = c('FDR', 'Significance')
  mylines = NULL
  
  i = 1
  for (l in labels) {
    info = infos[i, ]
    mylines = c(mylines, paste0( l,', ', 't(', ff(info['df']), ')=', ff(info['t value']), 
                                 ', p=', ff(info['Pr(>|t|)']),  ', p(FDR)=', ff(info['FDR']) ))
    
    i = i + 1
  }
  if (!output_text) mylines = NULL
  mylines = c(mylines, '---------------', "ROI,Estimate,Std. Error,df,t,p-value,Corrected p-value, Significance")
  
  i = 1
  #  for (g in groups) {
  for (l in labels) {
    info = infos[i, ]
    mylines = c(mylines, paste(l, ff(info['Estimate']), ff(info['Std. Error']), ff(info['df']), 
                               ff(info['t value']), ff(info['Pr(>|t|)']), ff(info['FDR']), info['Significance'], sep = ','))
    i = i + 1
  }
  #}
  writetofile(mylines, myfile)
  
}


do_tests_diff = function(mydata, par_diff, myfile, output_text = F){
  unlink(myfile)
  labels = sort(unique(mydata$label))
  groups = unique(mydata$GROUP)
  infos = NULL
  for (g in groups) {
    for (l in labels) {
      
      model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + (1 |subject), data = mydata %>% filter(label == l, GROUP == g))
      infos = rbind(infos, summary(model)$coefficients[par_diff, ])
      if(isSingular(model)) {
        print(myfile)
        print(summary(model))
      }
      
    }
    groupmodel = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + (1 |subject), data = mydata %>% filter( GROUP == g))
    info = summary(groupmodel)$coefficients[par_diff, ]
    myline = paste0(g, ', ', 't(', ff(info['df']), ')=', ff(info['t value']), ', p=', ff(info['Pr(>|t|)']) )
    writetofile(myline, myfile)
  }
  
  infos = cbind(infos, p.adjust(infos[, 'Pr(>|t|)'], method = 'fdr'))
  significance = ifelse(infos[, 5] < 0.05, '*', 'n. s.')
  significance[infos[, 6] < 0.05] = '**'
  infos = cbind(infos, significance)
  colnames(infos)[6:7] = c('FDR', 'Significance')
  mylines = NULL
  
  i = 1
  for (g in groups) {
    for (l in labels) {
      info = infos[i, ]
      mylines = c(mylines, paste0(g, ', ', l,', ', 't(', ff(info['df']), ')=', ff(info['t value']), 
                                  ', p=', ff(info['Pr(>|t|)']),  ', p(FDR)=', ff(info['FDR']) ))
      i = i + 1
    }
  }
  
  if (!output_text) mylines = NULL
  
  mylines = c(mylines, '---------------', "Group,ROI,Estimate,Std. Error,df,t,p-value,Corrected p-value, Significance")
  
  i = 1
  for (g in groups) {
    for (l in labels) {
      info = infos[i, ]
      mylines = c(mylines, paste(g, l, ff(info['Estimate']), ff(info['Std. Error']), ff(info['df']), 
                                 ff(info['t value']), ff(info['Pr(>|t|)']), ff(info['FDR']), info['Significance'], sep = ','))
      i = i + 1
    }
  }
  writetofile(mylines, myfile)
  
}


do_variability_analysis = function(meas, suffix0, suffix1, par = NULL, 
                                   par_diff = NULL, 
                                   ylimit = NULL, 
                                   ylimit_diff = NULL, 
                                   odd = "Trained Untrained",
                                   analysis_type = 'groupxtrainingxmeasure',
                                   remove_wave = NULL, 
                                   add_legend = F) {
  legend.position = c(0.12, 0.9)
  legend.position.diff = c(0.12, 0.1)
  suffix = paste0(suffix0, suffix1)
  data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_', suffix0, '.csv')
  print(paste(meas, suffix))
  #data_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_mask-cross.csv'
  #output_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/clf_acc.csv'
  output_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/', suffix, '.csv')
  figs_dir = '/data/lv0/MotorSkill/figs/variability'
  data = read.table(data_file, header = T, sep = ',') %>% arrange(subject, session, hemi, label)
  
  data$GROUP = "Intervention"
  data$GROUP[grep("sub-lue.2", data$subject)] = "Control"
  
  # remove particular  wave
  if (!is.null(remove_wave)) data = data[ -grep(paste0('sub-lue', remove_wave), data$subject), ]
  
  if (meas == 'alpha_PCM'){
    mymeasures = c("alpha_trained_PCM", "alpha_untrained_PCM")
    data$value = data$alpha_trained_PCM - data$alpha_untrained_PCM
  } 
  
  if (meas == 'alpha'){
    mymeasures = c("alpha_trained", "alpha_untrained")
    data$value = data$alpha_trained - data$alpha_untrained
    ylabel = "Variability index"
    ylabel_diff = "Difference in variability score"
    
    if (length(odd) == 1 & odd[1] == 'Trained Different') {
      formula = 'value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*MEASURE + (1 | subject)'
      contrast = c(0, 0, 0, 0, 0, 0, -1)
      remove_measure = ""
    }
    
    
  } 
  
  if (meas == 'valid'){
    mymeasures = c("valid") #, "clf_acc_trained_untrained")
    data$value = data$valid_runs
    ylabel = "Number of valid trials"
  } 
  
  if (meas == 'clf'){
    if (suffix1 == '-all') mymeasures = c("clf_acc_trained", "clf_acc_untrained", "clf_acc_trained_untrained")
    if (suffix1 == '-different') mymeasures = c("clf_acc_trained", "clf_acc_untrained")
    data$value = data$clf_acc_trained - data$clf_acc_untrained
    ylabel = "Classification accuracy"
    ylabel_diff = "Difference in classification accuracy"
    legend.position = c(0.12, 0.85)
    if (length(odd) == 1 & odd == c("Trained Untrained")) {
      formula = 'value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*MEASURE + (1 | subject)'
      contrast = c(0, 0, 0, 0, 0, 0, 1)
      remove_measure = "Trained"
    }
    
    if (length(odd) == 2) {
      formula = 'value ~ 1 + FD + SYSTEM + CONFIGURATION'
      contrast = c(1, 0, 0, 0)
      remove_measure = c("Trained", "Trained Untrained")
      data$clf_acc_untrained = data$clf_acc_untrained - 0.5
      data = subset(data, GROUP == 'Control')
    }
    
    clean_same = T
  } 
  
  if (meas == 'xeuclidean'){
    if (suffix1 == '-same') mymeasures = c("xeuclidean_trained_same", "xeuclidean_untrained_same")
    if (suffix1 == '-different') mymeasures = c("xeuclidean_trained_different", "xeuclidean_untrained_different")
    if (suffix1 == '-all') mymeasures = c("xeuclidean_trained_same", "xeuclidean_untrained_same", "xeuclidean_trained_different", "xeuclidean_untrained_different", "xeuclidean_trained_untrained")
    if (suffix1 == '-same') data$value = data$xeuclidean_trained_same - data$xeuclidean_untrained_same 
    else data$value = data$xeuclidean_trained_different - data$xeuclidean_untrained_different
    ylabel = "Squared euclidean distance"
    
  }
  
  if (meas == 'xnobis'){
    if (suffix1 == '-same') mymeasures = c("xnobis_trained_same", "xnobis_untrained_same")
    if (suffix1 == '-untrained') mymeasures = c("xnobis_untrained_same", "xnobis_untrained_different")
    if (suffix1 == '-trained') mymeasures = c("xnobis_trained_same", "xnobis_trained_different")
    if (suffix1 == '-different') mymeasures = c("xnobis_trained_different", "xnobis_untrained_different", "xnobis_trained_untrained")
    if (suffix1 == '-all') mymeasures = c("xnobis_trained_same", "xnobis_untrained_same", "xnobis_trained_different", "xnobis_untrained_different", "xnobis_trained_untrained")
    
    if (suffix1 == '-same') data$value = data$xnobis_trained_same - data$xnobis_untrained_same 
    if (suffix1 == '-trained') data$value = data$xnobis_trained_different - data$xnobis_trained_same
    if (suffix1 == '-untrained') data$value = data$xnobis_untrained_same - data$xnobis_untrained_different
    if (suffix1 == '-different') data$value = data$xnobis_trained_different - data$xnobis_untrained_different
    if (suffix1 == '-all') data$value = data$xnobis_trained_different - data$xnobis_untrained_different
    
    ylabel = "Cross-nobis dissimilarity"
    ylabel_diff = "Difference in cross-nobis dissimilarity\n(same - different type)"
    legend.position = 'bottom'
  } 
  
  if (meas == 'xnobisratio'){
    if (suffix1 == '-same') mymeasures = c("xnobisratio_trained_same", "xnobisratio_untrained_same")
    if (suffix1 == '-untrained') mymeasures = c("xnobisratio_untrained_same", "xnobisratio_untrained_different")
    if (suffix1 == '-trained') mymeasures = c("xnobisratio_trained_same", "xnobisratio_trained_different")
    if (suffix1 == '-different') mymeasures = c("xnobisratio_trained_different", "xnobisratio_untrained_different")
    if (suffix1 == '-all') mymeasures = c("xnobisratio_trained_same", "xnobisratio_untrained_same", "xnobisratio_trained_different", "xnobisratio_untrained_different", "xnobisratio_trained_untrained")
    
    if (suffix1 == '-same') data$value = data$xnobisratio_trained_same - data$xnobisratio_untrained_same 
    if (suffix1 == '-trained') data$value = data$xnobisratio_trained_different - data$xnobisratio_trained_same
    if (suffix1 == '-untrained') data$value = data$xnobisratio_untrained_different - data$xnobisratio_untrained_same
    if (suffix1 == '-different') data$value = data$xnobisratio_trained_different - data$xnobisratio_untrained_different
    if (suffix1 == '-all') data$value = data$xnobisratio_trained_different - data$xnobisratio_untrained_different
    
    ylabel = "Cross-nobis ratio"
  } 
  
  if (meas == 'xnobiscosine'){
    
    if (suffix1 == '-samediff') mymeasures = c("xnobiscosine_trained_same", "xnobiscosine_trained_different",
                                               "xnobiscosine_untrained_same", "xnobiscosine_untrained_different")
    if (suffix1 == '-same') mymeasures = c("xnobiscosine_trained_same", "xnobiscosine_untrained_same")
    if (suffix1 == '-untrained') mymeasures = c("xnobiscosine_untrained_same", "xnobiscosine_untrained_different")
    if (suffix1 == '-trained') mymeasures = c("xnobiscosine_trained_same", "xnobiscosine_trained_different")
    if (suffix1 == '-different') mymeasures = c("xnobiscosine_trained_different", "xnobiscosine_untrained_different")
    if (suffix1 == '-all') mymeasures = c("xnobiscosine_trained_same", "xnobiscosine_untrained_same", "xnobiscosine_trained_different", "xnobiscosine_untrained_different", "xnobiscosine_trained_untrained")
    
    if (suffix1 == '-samediff') data$value = data$xnobiscosine_same - data$xnobiscosine_different
    if (suffix1 == '-same') data$value = data$xnobiscosine_trained_same - data$xnobiscosine_untrained_same 
    if (suffix1 == '-trained') data$value = data$xnobiscosine_trained_different - data$xnobiscosine_trained_same
    if (suffix1 == '-untrained') data$value = data$xnobiscosine_untrained_different - data$xnobiscosine_untrained_same
    if (suffix1 == '-different') data$value = data$xnobiscosine_trained_different - data$xnobiscosine_untrained_different
    if (suffix1 == '-all') data$value = data$xnobiscosine_trained_different - data$xnobiscosine_untrained_different
    
    ylabel = "Cross-nobis cosine"
    
  } 
  
  if (meas == 'xnobisprod'){
    if (suffix1 == '-same') mymeasures = c("xnobisprod_trained_same", "xnobisprod_untrained_same")
    if (suffix1 == '-untrained') mymeasures = c("xnobisprod_untrained_same", "xnobisprod_untrained_different")
    if (suffix1 == '-trained') mymeasures = c("xnobisprod_trained_same", "xnobisprod_trained_different")
    if (suffix1 == '-different') mymeasures = c("xnobisprod_trained_different", "xnobisprod_untrained_different")
    if (suffix1 == '-all') mymeasures = c("xnobisprod_trained_same", "xnobisprod_untrained_same", "xnobisprod_trained_different", "xnobisprod_untrained_different", "xnobisprod_trained_untrained")
    
    if (suffix1 == '-same') data$value = data$xnobisprod_trained_same - data$xnobisprod_untrained_same 
    if (suffix1 == '-trained') data$value = data$xnobisprod_trained_different - data$xnobisprod_trained_same
    if (suffix1 == '-untrained') data$value = data$xnobisprod_untrained_different - data$xnobisprod_untrained_same
    if (suffix1 == '-different') data$value = data$xnobisprod_trained_different - data$xnobisprod_untrained_different
    if (suffix1 == '-all') data$value = data$xnobisprod_trained_different - data$xnobisprod_untrained_different
    
    ylabel = "Cross-nobis product"
  } 
  
  if (meas == 'xnobis_grouped'){
    
    if (length(odd) == 1 & odd[1] == 'Trained Different' & add_legend) {
      formula = 'value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*MEASURE + (1 | subject)'
      contrast = c(0, 0, 0, 0, 0, 0, 1)
      remove_measure = "Trained"
    }
    
    if (suffix1 == '-different') mymeasures = c("xnobis_grouped_trained_different", "xnobis_grouped_untrained_different", "xnobis_grouped_trained_untrained")
    if (suffix1 == '-different') data$value = data$xnobis_grouped_trained_different - data$xnobis_grouped_untrained_different
    
    ylabel = "Cross-nobis dissimilarity"
    ylabel_diff = "Difference in cross-nobis dissimilarity\n(same - different type)"
    legend.position = 'bottom'
    clean_same = T
    
  } 
  
  if (meas == 'xcosine'){
    if (suffix1 == '-same') mymeasures = c("xcosine_trained_same", "xcosine_untrained_same")
    if (suffix1 == '-different') mymeasures = c("xcosine_trained_different", "xcosine_untrained_different")
    if (suffix1 == '-all') mymeasures = c("xcosine_trained_same", "xcosine_untrained_same", "xcosine_trained_different", "xcosine_untrained_different", "xcosine_trained_untrained")
    if (suffix1 == '-same') data$value = data$xcosine_trained_same - data$xcosine_untrained_same 
    else data$value = data$xcosine_trained_different - data$xcosine_untrained_different
    ylabel = "Cosine (z-scored)"
  } 
  
  if (meas == 'xcorrelation'){
    if (suffix1 == '-same') mymeasures = c("xcorrelation_trained_same", "xcorrelation_untrained_same")
    if (suffix1 == '-untrained') mymeasures = c("xcorrelation_untrained_same", "xcorrelation_untrained_different")
    if (suffix1 == '-trained') mymeasures = c("xcorrelation_trained_same", "xcorrelation_trained_different")
    if (suffix1 == '-different') mymeasures = c("xcorrelation_trained_different", "xcorrelation_untrained_different",
                                                "xcorrelation_trained_untrained")
    if (suffix1 == '-all') mymeasures = c("xcorrelation_trained_same", "xcorrelation_untrained_same", "xcorrelation_trained_different", "xcorrelation_untrained_different", "xcorrelation_trained_untrained")
    
    if (suffix1 == '-same') data$value = data$xcorrelation_trained_same - data$xcorrelation_untrained_same 
    if (suffix1 == '-trained') data$value = data$xcorrelation_trained_same - data$xcorrelation_trained_different
    if (suffix1 == '-untrained') data$value = data$xcorrelation_untrained_same- data$xcorrelation_untrained_different
    if (suffix1 == '-different') data$value = data$xcorrelation_trained_different - data$xcorrelation_untrained_different
    if (suffix1 == '-all') data$value = data$xcorrelation_trained_different - data$xcorrelation_untrained_different
    
    #  data$value = data$correlation_trained_different - data$correlation_untrained_different 
    ylabel = "Correlation (z-scored)"
    ylabel_diff = "Difference in (z-scored) correlation\n(same - different type)"
    legend.position = 'bottom'
  } 
  
  if (meas == 'xcorrelation_unbiased'){
    if (suffix1 == '-different') mymeasures = c("xcorrelation_unbiased_trained_different", "xcorrelation_unbiased_untrained_different",
                                                "xcorrelation_unbiased_trained_untrained")
    if (suffix1 == '-different') data$value = data$xcorrelation_unbiased_trained_different - data$xcorrelation_unbiased_untrained_different
    
    ylabel = "Correlation (z-scored)"
    ylabel_diff = "Difference in (z-scored) correlation\n(trained - untrained type)"
    legend.position = 'bottom'
  } 
  
  if (meas == 'xcorrelation_grouped'){
    
    if (suffix1 == '-different') mymeasures = c("xcorrelation_grouped_trained_different", "xcorrelation_grouped_untrained_different",
                                                "xcorrelation_grouped_trained_untrained")
    if (suffix1 == '-different') data$value = data$xcorrelation_grouped_trained_different - data$xcorrelation_grouped_untrained_different
    
    ylabel = "Correlation (z-scored)"
    ylabel_diff = "Difference in (z-scored) correlation\n(trained - untrained)"
    legend.position = 'bottom'
  } 
  
  if (meas == 'xproduct_unbiased'){
    if (suffix1 == '-same') mymeasures = c("xproduct_unbiased_trained_same", "xproduct_unbiased_untrained_same")
    if (suffix1 == '-untrained') mymeasures = c("xproduct_unbiased_untrained_same", "xproduct_unbiased_untrained_different")
    if (suffix1 == '-trained') mymeasures = c("xproduct_unbiased_trained_same", "xproduct_unbiased_trained_different")
    if (suffix1 == '-different') mymeasures = c("xproduct_unbiased_trained_different", "xproduct_unbiased_untrained_different",
                                                "xproduct_unbiased_trained_untrained")
    if (suffix1 == '-all') mymeasures = c("xproduct_unbiased_trained_same", "xproduct_unbiased_untrained_same", "xproduct_unbiased_trained_different", "xproduct_unbiased_untrained_different", "xproduct_unbiased_trained_untrained")
    
    if (suffix1 == '-same') data$value = data$xproduct_unbiased_trained_same - data$xproduct_unbiased_untrained_same 
    if (suffix1 == '-trained') data$value = data$xproduct_unbiased_trained_same - data$xproduct_unbiased_trained_different
    if (suffix1 == '-untrained') data$value = data$xproduct_unbiased_untrained_same- data$xproduct_unbiased_untrained_different
    if (suffix1 == '-different') data$value = data$xproduct_unbiased_trained_different - data$xproduct_unbiased_untrained_different
    if (suffix1 == '-all') data$value = data$xproduct_unbiased_trained_different - data$xproduct_unbiased_untrained_different
    
    #  data$value = data$correlation_trained_different - data$correlation_untrained_different 
    ylabel = "Scalar product"
    ylabel_diff = "Difference in scalar product\n(trained - untrained)"
    legend.position = 'bottom'
  } 
  
  
  if (meas == 'xproduct_grouped'){
    if (suffix1 == '-same') mymeasures = c("xproduct_grouped_trained_same", "xproduct_grouped_untrained_same")
    if (suffix1 == '-untrained') mymeasures = c("xproduct_grouped_untrained_same", "xproduct_grouped_untrained_different")
    if (suffix1 == '-trained') mymeasures = c("xproduct_grouped_trained_same", "xproduct_grouped_trained_different")
    if (suffix1 == '-different') mymeasures = c("xproduct_grouped_trained_different", "xproduct_grouped_untrained_different",
                                                "xproduct_grouped_trained_untrained")
    if (suffix1 == '-all') mymeasures = c("xproduct_grouped_trained_same", "xproduct_grouped_untrained_same", "xproduct_grouped_trained_different", "xproduct_grouped_untrained_different", "xproduct_grouped_trained_untrained")
    
    if (suffix1 == '-same') data$value = data$xproduct_grouped_trained_same - data$xproduct_grouped_untrained_same 
    if (suffix1 == '-trained') data$value = data$xproduct_grouped_trained_same - data$xproduct_grouped_trained_different
    if (suffix1 == '-untrained') data$value = data$xproduct_grouped_untrained_same- data$xproduct_grouped_untrained_different
    if (suffix1 == '-different') data$value = data$xproduct_grouped_trained_different - data$xproduct_grouped_untrained_different
    if (suffix1 == '-all') data$value = data$xproduct_grouped_trained_different - data$xproduct_grouped_untrained_different
    
    #  data$value = data$correlation_trained_different - data$correlation_untrained_different 
    if (suffix1 == '-same') {
      if (add_legend) {
        formula = 'value ~ 1 + FD + SYSTEM + CONFIGURATION + GROUP*MEASURE + (1 | subject)'
        contrast = c(0, 0, 0, 0, 0, 0, -1)
        remove_measure = "Trained Untrained"
        
      }
      clean_same = T
      ylabel = "Cross-validated variance"
      ylabel_diff = "Difference in cross-validated variance estimate \n(trained - untrained)"
    } else {
      ylabel = "Scalar product"
      ylabel_diff = "Difference in scalar product"
      
    }
    legend.position = 'bottom'
    
  } 
  
  if (meas == 'mean_signal'){
    data$mean_signal_trained = data$mean_signal_trained/100
    data$mean_signal_untrained = data$mean_signal_untrained/100
    mymeasures = c("mean_signal_trained", "mean_signal_untrained")
    data$value = data$mean_signal_trained - data$mean_signal_untrained 
    ylabel = "% signal change"
    
  } 
  
  ##mymeasures = c("clf_acc_trained", "clf_acc_untrained")
  #c("mean_signal_trained", "mean_signal_untrained") #c("alpha_trained", "alpha_untrained")
  #mymeasures = c("crossnobis_trained", "crossnobis_untrained") #c("mean_signal_trained", "mean_signal_untrained") #c("alpha_trained", "alpha_untrained")
  
  # save it 
  #  output.data = data%>%
  #    group_by(label)%>%
  #    summarise(value=mean(value)) %>% arrange(-value) #, value_perm = mean(value_perm)
  
  #  write.table(output.data, file = output_file, sep = ',', col.names = T, row.names = F)
  
  nrois = 4
  # check data
  plot(sort(table(data$subject))/5/2, las = 2) # how many sessions per subject
  plot(sort(table(data$session))/5/2, las = 2) # how many sessions per timepoint
  
  # relabel rois and separate hemispheres
  control_labels = c("R_C1", "L_C1")
  data = data %>% filter( !label %in% control_labels)
  incomplete_subjects = c("sub-lue5207", "sub-lue3203") # too few timepoints
  data$label  = gsub("R_", "Right ", data$label)
  data$label  = gsub("L_", "Left ", data$label)
  data$label  = gsub("C1", "Control Region", data$label)
  data$uni_label  = gsub("Right ", "Average ", data$label)
  data$uni_label  = gsub("Left ", "Average ", data$uni_label)
  
  data = data %>% filter(!is.infinite(value))  %>% filter(!subject %in% incomplete_subjects, valid_runs >3) 
  
  data = data %>% group_by(subject) %>% 
    mutate(value=markoutliersIQR(value)) #%>% filter(!is.na(value)) %>% ungroup()
  
  data = merge(data, covars.table %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
               by = c("subject", "session", "GROUP"))
  
  data = left_join(data, motion %>% group_by(SUBJECT, TP) %>% 
                     summarise(FD = mean (FD)) %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
                   by = c("subject", "session"))
  #print(mymeasures)
  data = data %>% mutate(FD.valid = markoutliersIQR(FD)) %>% filter(!is.na(FD.valid))
  data.melt = reshape2::melt(data, 
                             id.vars = c("subject", "session", "GROUP", "hemi", "label", 
                                         "FD", "SYSTEM", "CONFIGURATION", "TRAINING", "TRAINING.Q", "TRAINING.A"),
                             variable.name = "MEASURE", 
                             value.name = "value") %>% mutate(SUBJECT = subject, value = as.numeric(value)) %>%
    filter(MEASURE %in% mymeasures) 
  
  
  # model the data + as.factor(CONFIGURATION)*as.factor(session)
  #
  
  # correction
  if (F) {
    data.melt$value.corr  = 0
    
    # for (l in unique(data.melt$label)) {
    #   print(l)
    #   model = lmer(clf_acc ~ 1 + FD + SYSTEM + (1 |subject), data = data %>% filter(label == l))
    #   print(summary(model))
    # }
    for (l in unique(data.melt$label)) {
      model = lmer(value ~ 1 + FD + SYSTEM + as.factor(CONFIGURATION) + GROUP*MEASURE*(TRAINING + TRAINING.Q) + 
                     + (1 + TRAINING + TRAINING.Q|SUBJECT), data = data.melt %>% filter(label == l))
      print(summary(model))
      data.melt$value.corr[data.melt$label == l] = predict(model, data.melt %>% filter(label == l) %>% 
                                                             mutate(FD = 0,
                                                                    SYSTEM = 'Classic',
                                                                    CONFIGURATION = '1'
                                                             )) + resid(model)
    }
    
  } else {
    data.melt$value.corr  = data.melt$value
    
  }
  
  
  # average
  data.mean = data%>%
    group_by(session, GROUP, label)%>%summarise(val.mean = mean(value, na.rm = T), 
                                                val.sem = sem(value)) 
  data.uni = data %>%
    group_by(subject, session, GROUP, uni_label)%>%summarise(value = mean(value, na.rm = T)) %>%
    group_by(session, GROUP, uni_label)%>%summarise(val.mean = mean(value, na.rm = T), 
                                                    val.sem = sem(value)) %>% mutate(label = uni_label)
  #data.mean = rbind(data.mean, data.uni)
  
  
  data.mean$ymin = data.mean$val.mean - data.mean$val.sem
  data.mean$ymax = data.mean$val.mean + data.mean$val.sem
  if (is.null(ylimit_diff)){
    ylimit_diff = c(min(data.mean$ymin), max(data.mean$ymax))
    ylimit_diff = ylimit_diff + diff(ylimit_diff)*c(-.05, .05)
  }
  myplot = ggplot(data.mean , aes(
    x = session - 1,
    y = val.mean,
    ymin = ymin, 
    ymax = ymax, 
    col = GROUP,
    group = GROUP
  )) + geom_line(size = 1) + 
    geom_point() + 
    geom_errorbar() + 
    facet_wrap(. ~ label, ncol = nrois) + 
    ylab(ylabel_diff) +
    xlab('Test session') + 
    ylim(ylimit_diff) + 
    geom_hline(yintercept = 0, size = 0.3, linetype = 2) + 
    scale_colour_manual(values = myPalette2) +
    scale_x_continuous(breaks = seq(0, 6)) +
                         theme_lh() + 
                         theme(legend.title = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               legend.text=element_text(size=14, face="bold"),
                               legend.box = 'vertical',
                               legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "pt"),
                               axis.text=element_text(size=18),
                               axis.title=element_text(size=20, face="bold"),
                               strip.text.x = element_text(size = 20, face="bold"),
                               strip.background = element_blank())
  
  
  if (add_legend) myplot = myplot + theme(legend.position = legend.position.diff)
  else myplot = myplot + theme(legend.position = "none")
  
  print(myplot)
  
  
  # hemi
  data.melt = data.melt%>% mutate(MEASURE = clean_measure_names(MEASURE, clean_same = clean_same))
  
  data.mean = data.melt%>%
    group_by(session, GROUP, label, MEASURE)%>%summarise(val.mean = mean(value.corr, na.rm = T), 
                                                         val.sem = sem(value.corr))%>% 
    mutate(group_measure = paste(GROUP, MEASURE)) 
  
  data.mean$ymin = data.mean$val.mean - data.mean$val.sem
  data.mean$ymax = data.mean$val.mean + data.mean$val.sem
  if (is.null(ylimit)){
    ylimit = c(min(data.mean$ymin), max(data.mean$ymax))
    ylimit = ylimit + diff(ylimit)*c(-.05, .05)
    
  }  
  
  if(length(table(data.mean$MEASURE)) == 2) myPalette = myPalette2
  myplot.all = ggplot(data.mean, aes(
    x = session - 1,
    y = val.mean,
    ymin = ymin, 
    ymax = ymax, 
    col = MEASURE,
    group = group_measure,
    linetype = GROUP
  )) + 
    geom_line(size = 1) + 
    geom_point() + 
    geom_errorbar() + 
    facet_wrap(. ~label, ncol = nrois) + 
    ylab(ylabel) +
    xlab('Test session') + 
    ylim(ylimit) + 
    theme_lh() + 
    scale_colour_manual(values = myPalette) +  
    scale_x_continuous(breaks= seq(0, 6)) +
    theme_lh() + theme(legend.position = c(0.8, 0.8), legend.title=element_blank()) + 
    theme(text = element_text(size = 9),
          strip.background = element_blank(),
          legend.title = element_blank(),
          legend.box = 'vertical',
          legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "pt"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text=element_text(size=16),
          axis.text=element_text(size=18),
          axis.title=element_text(size=20, face="bold"),
          strip.text.x = element_text(size = 20, face="bold"))
  
  
  
  if (!is.null(contrast)) {
    data.melt = data.melt %>% filter(! MEASURE %in% remove_measure)
    stats.timepoints = test_model(data.melt, seq(7), 
                                  formula, meas, 
                                  contrast = contrast, print_model = F) %>% mutate(label = Label)
    
    #print(stats.timepoints)
    View(stats.timepoints)
    stats.timepoints.file = file.path(figs_dir, paste(meas, suffix, 'stats_timepoints.txt', sep = '-'))
    write.table(stats.timepoints[seq(11)], file = stats.timepoints.file, sep = ',', col.names = T, row.names = F, quote = F)
    
    myplot.all = myplot.all + geom_text( data = stats.timepoints, aes(x = Session - 1, 
                                                                      y = ylimit[2]*0.95, label = Significance), inherit.aes = F, size = STARSIZE) 
    
  }  
  
  if (add_legend) myplot.all = myplot.all + theme(legend.position = legend.position)
  else myplot.all = myplot.all + theme(legend.position = "none")
  
  print(myplot.all)
  
  ggsave(
    file.path(figs_dir, paste(meas, suffix, 'difference.png', sep = '-')),
    plot = myplot,
    width = WIDTH,
    height = HEIGHT,
    dpi = DPI,
    units = 'cm'
  )
  
  ggsave(
    file.path(figs_dir, paste(meas, suffix, 'all.png', sep = '-')),
    plot = myplot.all,
    width = WIDTH,
    height = HEIGHT,
    dpi = DPI,
    units = 'cm'
  )
  if (!is.null(par)) do_tests(data.melt, par, file.path(figs_dir, paste(meas, suffix, gsub(" ", "", par), 'tests.txt', sep = '-')), odd = odd,
                              analysis_type = analysis_type)
  
  if (!is.null(par_diff)) do_tests_diff(data, par_diff, file.path(figs_dir, paste(meas, suffix, 'tests_diff.txt', sep = '-'))) 
  
}

test_model = function(data, test_sessions, formula, name, contrast, print_model = F, method = "fdr", family = NULL){
  stats = NULL
  singular = NULL
  hasconverged = NULL
  #  browser()
  data$MEASURE = relevel(factor(data$MEASURE), ref = 'Untrained')
  #  strict_tol <- lmerControl(optCtrl=list(method="nlminb"), optimizer = "optimx")
  labels = unique(data$label)
  for (mylabel in labels){
    
    for (sess in test_sessions){
      
      mydata = subset(data, session == sess & label == mylabel)
      if (length(unique(mydata$SYSTEM)) == 1)
      {
        contrast.1 = contrast[-3]
        formula.1 = gsub('SYSTEM \\+ ', '', formula)
      } else {
        contrast.1 = contrast
        formula.1 = formula
      }
      if(length(grep("\\|", formula))==1) model = lmer( formula.1, data = mydata)
      else model = lm( formula.1, data = mydata)
      #      browser()
      #print(
      #  ggplot(subset(data, sess_num == sess), aes(x = seq_train, y = log(meanMT), col = seq_train, group = seq_train)) + geom_boxplot() + facet_grid(CONFIGURATION.SIMPLE ~  group)
      #)
      if (print_model) print(summary(model))
      if(length(grep("\\|", formula))==1) {
        singular = c(singular, isSingular(model))
        hasconverged = c(hasconverged, model@optinfo$conv$opt)
      } else {
        singular = F
        hasconverged = T
      }
      g = glht(model, linfct = rbind(contrast.1), alternative = "greater")
      gg = summary(model)$coefficients[abs(contrast.1) == 1, ]
      pvalue = unlist(summary(g)$test[ 'pvalues'])
      
      if (!is.null(family)){
        gg = c(gg[1:2], NA, gg[3:4])
      }
      gg[5] = pvalue
      stats = rbind(stats, c(mylabel, sprintf("%d#", sess-1), sprintf("%d#", (sess-1)*6), gg))  
    }
    print("==========================")
    
    if (nlevels(as.factor(data$GROUP)) == 1){
      model.ROI = lmer('value ~ 1 + FD + SYSTEM + CONFIGURATION + (1|subject)', data = subset(data, label == mylabel))
      print(summary(model.ROI))
      print("==========================")
    }
    
  }
  stats = as.data.frame(stats)
  colnames(stats) = c('Label', 'MRI session', 'Behavioral session', 'Estimate', 'Std. Error', 'df', 't', 'p(uncorrected)')
  
  stats$`p(corrected)` = p.adjust(stats$p, method = method)
  stats = stats %>% group_by(Label) %>% mutate(`p(corrected_label)` = p.adjust(`p(uncorrected)`, method = method))
  stats.format = as.data.frame(apply(stats, c(1,2), ff))
  stats.format$Significance = ifelse(stats$`p(uncorrected)` < 0.05, '*', ' ')
  #  stats.format$Significance[ stats$`p(uncorrected)` < 0.05 ] =  '*'
  stats.format$Significance[ stats$`p(corrected_label)` < 0.05 ] = '**'
  stats.format$Significance[ stats$`p(corrected)` < 0.05 ] = '***'
  stats.format$Name = name
  stats.format$Session = test_sessions  
  stats.format$Singular = singular  
  stats.format$converged = hasconverged  
  return(stats.format)
}


if (F) {
  WIDTH = 30; HEIGHT = 24; DPI = 1000
  meas = 'xcorrelation'
  #meas = 'xcosine'
  #meas = 'xeuclidean'
  #  meas = 'xnobisratio'
  #  meas = 'xnobisprod'
  #  meas = 'xnobis'
  #  meas = 'xnobiscosine'
  #meas = 'mean_signal'
  #meas = 'clf'
  #meas = 'valid'
  #suffix0 = 'mask-cross-perm'
  #suffix0 = 'mask-cross-derivatives'
  suffix0 = 'mask-cross-runprew'
  
  #  suffix0 = 'mask-cross'
  #suffix0 = 'mask-cross-runprew-perm'
  #suffix0 = 'mask-cross-noprew'
  suffix1 = '-all'
  #suffix1 = '-same'
  suffix1 = '-samediff'
  suffix1 = '-different'
  #  suffix1 = '-untrained'
  #suffix1 = '-trained'
  do_variability_analysis(meas, suffix0, suffix1)
}
