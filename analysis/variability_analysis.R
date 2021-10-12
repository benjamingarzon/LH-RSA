library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(lme4)
library(lmerTest)

# remove motion outliers...?
# correct for FD. FD vs type of 
# correct for number of correct trials

do_tests = function(mydata, par, myfile, odd, analysis_type = 'groupxtrainingxmeasure', output_text = F){
  
  unlink(myfile)
  mydata = mydata %>% filter(! MEASURE %in% odd)
  labels = sort(unique(mydata$label))
  groups = unique(mydata$GROUP)
  infos = NULL
  mylines = c(par)
  for (l in labels) {
    if (analysis_type == 'groupxtrainingxmeasure') 
      model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + (GROUP*TRAINING*MEASURE) + (1 |subject), data = mydata %>% 
                     filter(label == l))
    if (analysis_type == 'groupxtraining') 
        model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + (GROUP*TRAINING) + (1 |subject), data = mydata %>% 
                     filter(label == l))
    if (analysis_type == 'training') 
      model = lmer(value ~ 1 + FD + SYSTEM + CONFIGURATION + TRAINING + (1 |subject), data = mydata %>% 
                     filter(label == l))
    if(isSingular(model)) {
      print(myfile)
      print(summary(model))
    }
    infos = rbind(infos, summary(model)$coefficients[par, ])
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
                                   remove_wave = NULL) {
  suffix = paste0(suffix0, suffix1)
  data_file = paste0('/data/lv0/MotorSkill/fmriprep/analysis/surf/roi_scores_', suffix0, '.csv')
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
    ylabel = "Variability score"
  } 
  
  if (meas == 'valid'){
    mymeasures = c("valid") #, "clf_acc_trained_untrained")
    data$value = data$valid_runs
    ylabel = "Number of valid trials"
  } 
  
  if (meas == 'clf'){
    mymeasures = c("clf_acc_trained", "clf_acc_untrained") #, "clf_acc_trained_untrained")
    data$value = data$clf_acc_trained - data$clf_acc_untrained
    ylabel = "Classification accuracy"
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
    ylabel_diff = "Difference in cross-nobis dissimilarity (same - different type)"
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
    ylabel_diff = "Difference in (z-scored) correlation  (same - different type)"
    
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
  incomplete_subjects = c("sub-lue5207", "sub-lue3203")
  data$label  = gsub("R_", "Right ", data$label)
  data$label  = gsub("L_", "Left ", data$label)
  data$label  = gsub("C1", "Control Region", data$label)
  data$uni_label  = gsub("Right ", "Average ", data$label)
  data$uni_label  = gsub("Left ", "Average ", data$uni_label)
  
  data = data %>% filter(!is.infinite(value))  %>% filter(!subject %in% incomplete_subjects, valid_runs >3) 
  
  data = data %>% group_by(subject) %>% 
    mutate(value=markoutliersIQR(value)) %>% 
    filter(!is.na(value)) %>% ungroup()
  
  data = merge(data, covars.table %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
               by = c("subject", "session", "GROUP"))
  
  data = left_join(data, motion %>% group_by(SUBJECT, TP) %>% 
                     summarise(FD = mean (FD)) %>% mutate(subject = paste0('sub-', SUBJECT), session = TP), 
                   by = c("subject", "session"))
  #print(mymeasures)
  
  data.melt = reshape2::melt(data, 
                             id.vars = c("subject", "session", "GROUP", "hemi", "label", 
                                         "FD", "SYSTEM", "CONFIGURATION", "TRAINING", "TRAINING.Q", "TRAINING.A"),
                             variable.name = "MEASURE", 
                             value.name = "value") %>% mutate(SUBJECT = subject, value = as.numeric(value)) %>%
    filter(MEASURE %in% mymeasures) 
  
  # relabel measures 
  
  
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
      print(l)
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
    x = session,
    y = val.mean,
    ymin = ymin, 
    ymax = ymax, 
    col = GROUP,
    group = GROUP
  )) + geom_line() + 
    geom_point() + 
    geom_errorbar() + 
    facet_wrap(. ~ label, ncol = nrois) + 
    ylab(ylabel_diff) +
    xlab('Test session') + 
    ylim(ylimit_diff) + 
    geom_hline(yintercept = 0, size = 0.3, linetype = 2) + 
    theme_lh() + 
    scale_colour_manual(values = myPalette2) +
    theme(legend.title = element_blank(), legend.position = 'bottom',  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text=element_text(size=18),
          axis.text=element_text(size=18),
          axis.title=element_text(size=20, face="bold"),
          strip.text.x = element_text(size = 20, face="bold"))
  
  
  
  print(myplot)
  
  
  # hemi
  data.melt = data.melt%>% mutate(MEASURE = clean_measure_names(MEASURE))
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
    x = session,
    y = val.mean,
    ymin = ymin, 
    ymax = ymax, 
    col = MEASURE,
    group = group_measure,
    linetype = GROUP
  )) + 
    geom_line(size = 0.6) + 
    geom_point() + 
    geom_errorbar() + 
    facet_wrap(. ~label, ncol = nrois) + 
    ylab(ylabel) +
    xlab('Test session') + 
    ylim(ylimit) + 
    theme_lh() + 
    scale_colour_manual(values = myPalette) + 
    theme(text = element_text(size = 9),
          legend.title = element_blank(), legend.position = 'bottom', 
          legend.box = 'horizontal',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text=element_text(size=18),
          axis.text=element_text(size=18),
          axis.title=element_text(size=20, face="bold"),
          strip.text.x = element_text(size = 20, face="bold"))
  
  
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
