
##### remove zero cols from computations?

rm(list = ls())

setwd("~/Software/LeftHand/analysis")
source('./plot_funcs.R')
source('./load_covariates.R')
source('./variability_analysis.R')
#source('./session_variability_analysis.R')

########################
# Figure options
########################
WIDTH = 30; HEIGHT = 24; DPI = 1000

########################
# Variability analysis
########################
# some selected figures
for (suffix0 in c('mask-cross-runprew', 'mask-cross-runprew-perm')) { #
  if (F){
 #do_variability_analysis('clf', suffix0, '-all', ylimit = c(0.45, 0.70), add_legend = suffix0 == 'mask-cross-runprew-perm')    
 do_variability_analysis('clf', suffix0, '-all', ylimit = c(0.45, 0.70), 
                         odd = c('Trained Different'))    
 do_variability_analysis('clf', suffix0, '-all', ylimit = c(0.45, 0.70), 
                         odd = c('Trained Different', 'Trained Untrained'))    
 do_variability_analysis('clf', suffix0, '-all', ylimit = c(0.45, 0.70), 
                         add_legend = suffix0 == 'mask-cross-runprew-perm')
 
 do_variability_analysis('alpha', suffix0, '-all', odd = c('Trained Different'), add_legend = T)     
# do_variability_analysis('alpha', suffix0, '-all', par = c('MEASUREUntrained'), add_legend = T)  
  }
  do_variability_analysis('xnobis_grouped', suffix0, '-different', par = c('TRAINING'),
                          odd = c('Trained Different', 'Untrained Different'),
                          analysis_type = 'groupxtraining')

  do_variability_analysis('xnobis_grouped', suffix0, '-different', par = c('GROUPIntervention'),
                          odd = c('Trained Different', 'Untrained Different'),
                          analysis_type = 'groupxtraining')

  do_variability_analysis('xnobis_grouped', suffix0, '-different', par = c('GROUPIntervention:TRAINING'),
                          odd = c('Trained Different', 'Untrained Different'),
                          analysis_type = 'groupxtraining') #xproduct_grouped

  do_variability_analysis('xnobis_grouped', suffix0, '-different', par = c('GROUPIntervention:TRAINING'),
                         odd = c('Trained Different'),
                         analysis_type = 'groupxtraining', add_legend = T) #xproduct_grouped
  
  do_variability_analysis('xproduct_grouped', suffix0, '-same', par = c('GROUPIntervention'), analysis_type = 'groupxmeasure')
  do_variability_analysis('xproduct_grouped', suffix0, '-same', par = c('MEASUREUntrained'), analysis_type = 'groupxmeasure')
  do_variability_analysis('xproduct_grouped', suffix0, '-same', par = c('GROUPIntervention:MEASUREUntrained'),
                          add_legend = T, analysis_type = 'groupxmeasure')
  

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
# stophere
# ########################
# # Session variability analysis
# ########################
# mylabels = c('Right PM', 'Right SPL', 'Right PS') #'Right Control Region', , 'Right SMA'
# #, 'Left PM', 'Left Control Region',  'Left SPL', 'Left PS', 'Left SMA')
# 
# for (suffix0 in c('mask-cross', 'mask-cross-perm')) { #mask-cross-derivatives-perm, 'mask-cross-perm', , 'mask-cross-derivatives'
#   for (measure in c('xcosine', 'xcorrelation'))  { #'xnobis', 
#     do_session_variability_analysis(measure, suffix0, '-all', mylabels)
#     do_session_variability_analysis(measure, suffix0, '-same', mylabels)
#     do_session_variability_analysis(measure, suffix0, '-different', mylabels)
#   }
# }
# 
#   
# 
# # do_variability_analysis('clf', suffix0, '-all', ylimit = c(0.45, 0.70), add_legend = suffix0 == 'mask-cross-runprew-perm')    
# # do_variability_analysis('clf', suffix0, '-different', ylimit = c(0.45, 0.65))    
# # do_variability_analysis('alpha', suffix0, '-all', par = c('GROUPIntervention'), add_legend = T, analysis_type = 'groupxtraining')  
# # do_variability_analysis('alpha', suffix0, '-all', par = c('MEASUREUntrained'), add_legend = T)  
# 
# # do_variability_analysis('xcorrelation_grouped', suffix0, '-different', add_legend = T)
# do_variability_analysis('xnobis_grouped', suffix0, '-different', add_legend = T)  
# #  do_variability_analysis('xproduct_grouped', suffix0, '-same', add_legend = T)  
# # do_variability_analysis('xcorrelation_grouped', suffix0, '-different', add_legend = T)  
# 
# #  do_variability_analysis('xproduct_unbiased', suffix0, '-different', add_legend = T)  
# #do_variability_analysis('xcorrelation_unbiased', suffix0, '-different', add_legend = T)  
# 
# 
# #do_variability_analysis('alpha', suffix0, '-all', par = c('GROUPIntervention:MEASUREUntrained'), add_legend = T, analysis_type = 'groupxmeasure')  
# #do_variability_analysis('alpha', suffix0, '-all', par = c('MEASUREUntrained'), add_legend = T, analysis_type = 'groupxmeasure')  
# 
# if (T){
#   do_variability_analysis('xcorrelation', suffix0, '-same', par = c('GROUPIntervention'), analysis_type = 'groupxmeasure')
#   do_variability_analysis('xcorrelation', suffix0, '-same', par = c('MEASUREUntrained Same'), analysis_type = 'groupxmeasure')
#   do_variability_analysis('xcorrelation', suffix0, '-same', par = c('GROUPIntervention:MEASUREUntrained Same'), ylimit = c(0.05, 0.13), add_legend = T, analysis_type = 'groupxmeasure')
#   do_variability_analysis('xcorrelation', suffix0, '-different', par = c('GROUPIntervention'), analysis_type = 'groupxmeasure')
#   do_variability_analysis('xcorrelation', suffix0, '-different', par = c('MEASUREUntrained Different'), analysis_type = 'groupxmeasure')
#   do_variability_analysis('xcorrelation', suffix0, '-different', par = c('GROUPIntervention:MEASUREUntrained Different'), ylimit = c(0.05, 0.13), add_legend = T, analysis_type = 'groupxmeasure')
#   
#   do_variability_analysis('xcorrelation', suffix0, '-untrained', par_diff = '(Intercept)', 
#                           ylimit_diff = c(-2e-3, 5e-3), add_legend = T)
#   
#   do_variability_analysis('xnobis', suffix0, '-same', 
#                           par = c('GROUPIntervention'),
#                           odd = c('Trained Same'), 
#                           analysis_type = 'groupxtraining')
#   
#   do_variability_analysis('xnobis', suffix0, '-same', 
#                           par = c('TRAINING'),
#                           odd = c('Trained Same'), 
#                           analysis_type = 'training')
#   
#   do_variability_analysis('xnobis', suffix0, '-same', 
#                           par = c('GROUPIntervention:TRAINING'), 
#                           odd = c('Trained Same'), 
#                           analysis_type = 'groupxtraining',
#                           ylimit = c(-0.01, 0.13), add_legend = T)
#   
#   do_variability_analysis('xnobis', suffix0, '-different', par = c('TRAINING'), 
#                           odd = c('Trained Different', 'Untrained Different'), 
#                           analysis_type = 'training')
#   
#   do_variability_analysis('xnobis', suffix0, '-different', par = c('GROUPIntervention'), 
#                           odd = c('Trained Different', 'Untrained Different'),
#                           analysis_type = 'groupxtraining', add_legend = T)
#   
#   do_variability_analysis('xnobis', suffix0, '-different', par = c('GROUPIntervention:TRAINING'), 
#                           odd = c('Trained Different', 'Untrained Different'),
#                           analysis_type = 'groupxtraining',
#                           ylimit = c(-0.01, 0.13), add_legend = T)
#   
#   do_variability_analysis('xnobis', suffix0, '-different', par = c('GROUPIntervention'), 
#                           odd = c('Trained Different', 'Untrained Different'),
#                           analysis_type = 'groupxtraining', add_legend = T)
#   
#   do_variability_analysis('xnobis', suffix0, '-untrained', par_diff = '(Intercept)', 
#                           ylimit_diff = c(-.08, .04), add_legend = T)
# }
# 
# }
# 
# 
# 
# stophere
# if(F) {
#   
#   for (suffix0 in c('mask-cross', 'mask-cross-perm')) { #mask-cross-derivatives-perm, 'mask-cross-perm',
#     print(suffix0)
#     for (measure in c( 'xnobis', 'xcosine', 'xcorrelation', 'xeuclidean', 'mean_signal', 'alpha', 'clf')){
#       #  for (measure in c( 'clf')){
#       do_variability_analysis(measure, suffix0, '-all')
#       
#     }
#     for (measure in c('xnobis', 'cosine', 'correlation'))  {
#       do_variability_analysis(measure, suffix0, '-different')
#       do_variability_analysis(measure, suffix0, '-same')
#     }
#     
#   }
#   
# }
