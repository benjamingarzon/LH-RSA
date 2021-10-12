
##### remove zero cols from computations?

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
# some selected figures
for (suffix0 in c('mask-cross-runprew', 'mask-cross-runprew-perm')) {
  if (T){
 do_variability_analysis('xcorrelation', suffix0, '-same', par = c('GROUPIntervention'), analysis_type = 'groupxtraining')
 do_variability_analysis('xcorrelation', suffix0, '-same', par = c('MEASUREUntrained Same'))
 do_variability_analysis('xcorrelation', suffix0, '-same', par = c('GROUPIntervention:MEASUREUntrained Same'))
 do_variability_analysis('xcorrelation', suffix0, '-different', par = c('GROUPIntervention'), analysis_type = 'groupxtraining')
 do_variability_analysis('xcorrelation', suffix0, '-different', par = c('MEASUREUntrained Different'))
 do_variability_analysis('xcorrelation', suffix0, '-different', par = c('GROUPIntervention:MEASUREUntrained Different'))
  }
  
 do_variability_analysis('xcorrelation', suffix0, '-untrained', par_diff = '(Intercept)', 
                          ylimit_diff = c(-2e-3, 5e-3))
  
 do_variability_analysis('xnobis', suffix0, '-same', 
                         par = c('TRAINING'),
                         odd = c('Trained Same'), 
                         analysis_type = 'training')

 do_variability_analysis('xnobis', suffix0, '-same', 
                        par = c('GROUPIntervention:TRAINING'), 
                        odd = c('Trained Same'), 
                        analysis_type = 'groupxtraining')

  do_variability_analysis('xnobis', suffix0, '-different', par = c('TRAINING'), 
                          odd = c('Trained Different', 'Untrained Different'), 
                          analysis_type = 'training')
  do_variability_analysis('xnobis', suffix0, '-different', par = c('GROUPIntervention:TRAINING'), 
                          odd = c('Trained Different', 'Untrained Different'),
                          analysis_type = 'groupxtraining')
#  do_variability_analysis('xnobis', suffix0, '-different', par = c('GROUPIntervention'), 
#                          odd = c('Trained Different', 'Untrained Different'),
#                          analysis_type = 'groupxtraining'))

  do_variability_analysis('xnobis', suffix0, '-untrained', par_diff = '(Intercept)', 
                          ylimit_diff = c(-.08, .04))
  
#  do_variability_analysis('xnobis', suffix0, '-different', par = c('MEASUREUntrained Different'), odd = 'Trained Different')
#  do_variability_analysis('xnobis', suffix0, '-different', par = c('GROUPIntervention:MEASUREUntrained Different'), odd = 'Trained Different')
#  do_variability_analysis('xnobis', suffix0, '-different', par = c('MEASURETrained Untrained'), odd = 'Untrained Different')
#  do_variability_analysis('xnobis', suffix0, '-different', par = c('GROUPIntervention:MEASURETrained Untrained'), odd = 'Untrained Different')
}


stophere
if(F) {

for (suffix0 in c('mask-cross', 'mask-cross-perm')) { #mask-cross-derivatives-perm, 'mask-cross-perm',
  print(suffix0)
  for (measure in c( 'xnobis', 'xcosine', 'xcorrelation', 'xeuclidean', 'mean_signal', 'alpha', 'clf')){
#  for (measure in c( 'clf')){
      do_variability_analysis(measure, suffix0, '-all')
    
  }
  for (measure in c('xnobis', 'cosine', 'correlation'))  {
    do_variability_analysis(measure, suffix0, '-different')
    do_variability_analysis(measure, suffix0, '-same')
  }
  
}

}
stophere
########################
# Session variability analysis
########################
mylabels = c('Right PM', 'Right SPL', 'Right PS') #'Right Control Region', , 'Right SMA'
#, 'Left PM', 'Left Control Region',  'Left SPL', 'Left PS', 'Left SMA')

for (suffix0 in c('mask-cross', 'mask-cross-perm')) { #mask-cross-derivatives-perm, 'mask-cross-perm', , 'mask-cross-derivatives'
  for (measure in c('xcosine', 'xcorrelation'))  { #'xnobis', 
    do_session_variability_analysis(measure, suffix0, '-all', mylabels)
    do_session_variability_analysis(measure, suffix0, '-same', mylabels)
    do_session_variability_analysis(measure, suffix0, '-different', mylabels)
  }
}

  

