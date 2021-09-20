rm(list = ls())

setwd("~/Software/LeftHand/analysis")
source('./plot_funcs.R')
source('./load_covariates.R')
source('./variability_analysis.R')
source('./session_variability_analysis.R')

########################
# Figure options
########################
WIDTH = 30; HEIGHT = 24; DPI = 1000

########################
# Variability analysis
########################
if(T) {

for (suffix0 in c('mask-cross', 'mask-cross-derivatives')) { #mask-cross-derivatives-perm, 'mask-cross-perm',
  print(suffix0)
#  for (measure in c( 'xnobis', 'cosine', 'correlation', 'euclidean', 'mean_signal', 'alpha', 'clf')){
  for (measure in c( 'clf')){
      do_variability_analysis(measure, suffix0, '-all')
    
  }
  for (measure in c('xnobis', 'cosine', 'correlation'))  {
    do_variability_analysis(measure, suffix0, '-different')
    do_variability_analysis(measure, suffix0, '-same')
  }
  
}

}

########################
# Session variability analysis
########################
mylabels = c('Right PM', 'Right SPL', 'Right PS') #'Right Control Region', , 'Right SMA'
#, 'Left PM', 'Left Control Region',  'Left SPL', 'Left PS', 'Left SMA')

for (suffix0 in c('mask-cross')) { #mask-cross-derivatives-perm, 'mask-cross-perm', , 'mask-cross-derivatives'
  for (measure in c('xcosine', 'xcorrelation'))  { #'xnobis', 
    do_session_variability_analysis(measure, suffix0, '-all', mylabels)
    do_session_variability_analysis(measure, suffix0, '-same', mylabels)
    do_session_variability_analysis(measure, suffix0, '-different', mylabels)
  }
}

  

