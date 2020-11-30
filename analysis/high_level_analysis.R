rm(list = ls())
library(pracma)
library(ggplot2)
library(dplyr)

source('~/Software/ImageLMMR/ImageLMMR.R')
setwd("~/Software/LeftHand/analysis")
source('./mytests.R')
source('./load_covariates.R')

nifti_convert = '~/Software/workbench/bin_rh_linux64/wb_command -metric-convert -to-nifti'

collect_data = T

# take run into account
# define vars
# covariates
# check contrasts
# add variability contrast 
# covariates, coil type


WD = '~/Data/LeftHand/Lund1/fmriprep/analysis'

#schedules = read.table('~/Data/LeftHand/Lund1/responses/all_schedule_tables.csv', header = T)
rh.gii = '~/Data/LeftHand/Lund1/fmriprep/freesurfer/fsaverage6/surf/rh.midthickness.surf.gii' 
lh.gii = '~/Data/LeftHand/Lund1/fmriprep/freesurfer/fsaverage6/surf/lh.midthickness.surf.gii'

#contrastname = "UntrainedCorrect"
#contrastnum = 2

contrastname = "UntrainedCorrect_TrainedCorrect"
contrastnum = 4

#contrastname = "Spread"
analysis_type = 'surfR'  #volume, surfR/L spreadR/L  

NPROCS = 30
# list files, adapt depending on type of analysis
contrast = switch(which(analysis_type == c('volume', 'surfR', 'surfL', 'spreadR', 'spreadL')),
                  paste0("/cope", contrastnum, ".nii.gz"),
                  paste0("/cope", contrastnum, ".func.gii"),
                  paste0("/cope", contrastnum, ".func.gii"),#                  paste0("\\bcope", contrastnum, "\\.func\\.gii$"),

                  #                  "\\brh\\.sl_within_spread_correlation_15\\.0\\.func\\.gii$",
#                  "\\blh\\.sl_within_spread_correlation_15\\.0\\.func\\.gii$"
                  "\\brh\\.sl_acc_svm_15\\.0\\.func\\.gii$",
                  "\\blh\\.sl_acc_svm_15\\.0\\.func\\.gii$"
)


if (analysis_type == 'volume') mysurf = '' else mysurf = ifelse(grepl("R", analysis_type, fixed = T), lh.gii, rh.gii)

image.list = system(sprintf("find %s/* | grep %s", WD, contrast), intern = TRUE) 
if (analysis_type %in% c('volume', 'surfR', 'surfL')) image.list =  image.list[grep(analysis_type, image.list)]
    
print(image.list)

# output files for analysis
DATADIR = file.path(WD, "higherlevel", contrastname)
dir.create(DATADIR)
dir.create(file.path(DATADIR, analysis_type))
dir.create(file.path(DATADIR, analysis_type, "tests"))

image_list_file = file.path(analysis_type, "image_list.txt")
write.table(image.list, file = file.path(DATADIR, image_list_file), row.names = F, quote = F, col.names = F)
setwd(DATADIR)


if (collect_data) {
  if (analysis_type == 'volume') {
   print("Converting and merging files")
   setwd(WD)
   system("fslmerge -t all_masks `find .. | grep bold_space-MNI152NLin2009cAsym_brainmask`")
   system("fslmaths all_masks -Tmean -mul 100 higherlevel/mask_sum") #`fslnvols      higherlevel/mask_mean all_masks`
   system("fslmaths all_masks -Tmin higherlevel/mask_min")
   system("mri_vol2vol --mov /home/benjamin.garzon/Data/LeftHand/Lund1/cat12/mask.nii.gz --targ higherlevel/mask_min.nii.gz --out higherlevel/GMprob.nii.gz --regheader")
   system(sprintf("fslmaths higherlevel/GMprob -thr 0.2 -bin -mul higherlevel/mask_min higherlevel/mask_min; mv higherlevel/mask_min.nii.gz higherlevel/%s/volume/", contrastname))
   system(sprintf("fslmaths higherlevel/GMprob -thr 0.2 -bin -mul higherlevel/mask_sum higherlevel/mask_sum; mv higherlevel/mask_sum.nii.gz higherlevel/%s/volume/", contrastname))
   system("rm all_masks.nii.gz higherlevel/GMprob.nii.gz")
  
   setwd(file.path(DATADIR, "volume"))
   system("rm mask.nii.gz; ln -s mask_min.nii.gz mask.nii.gz")
   system("fslmerge -t images `cat image_list.txt`")
   system("flirt -in mask -ref mask -out mask_4mm -interp nearestneighbour -applyisoxfm 4")
   system("flirt -in images -ref mask_4mm -out images_4mm -interp nearestneighbour -applyisoxfm 4")
   
  } else {
    setwd(DATADIR)  
    print("Converting and merging files")
    for (filename in image.list) {
      #print(filename)
      system(paste(nifti_convert, 
      filename, 
      paste0(filename, ".nii.gz")
      ))
      
      if (! analysis_type %in% c('volume', 'surfR', 'surfL')) system(sprintf("fslroi %s %s 0 1", filename, filename))
    }
    
    
    system(sprintf("IM=`cat %s | sed 's/.gii/.gii.nii.gz/g'`; fslmerge -t %s/images $IM; rm $IM", 
                   image_list_file, analysis_type))
    
    system(sprintf("fslmaths %s/images -Tstd -bin %s/mask", 
                   analysis_type, analysis_type))
  }
}

extract_substr = function(x, pattern, start, stop){
  i = strfind(x, pattern)
  if (!is.null(i)) substr = substr(x, i + start, i + stop)
  else substr = "-1"
  return(substr)
  
} 

doit = function(WD, IMAGES, MYTEST, OD, 
                MASK = '', 
                IMAGES_NAME = '',
                IMAGING_NAME = '',
                to_gifti = '', flip = F){
  
  setwd(WD)
  IMAGING_FILE = file.path(WD, IMAGING_NAME)
  MASK_FILE = file.path(WD, MASK)
  OUTPUT_DIR = file.path(WD, OD)
  # create regressors
  subject = sapply(IMAGES, extract_substr, "sub-", 4, 10)
  sess_num = sapply(IMAGES, extract_substr, "ses-", 4, 4) 
  run = sapply(IMAGES, extract_substr, "run", 3, 3) 
  
  DATA = data.frame(NAME = unlist(IMAGES), SUBJECT = subject, TP = sess_num, RUN = run) 
  DATA$SUBJECT.NUM = as.numeric(DATA$SUBJECT)
  
  #design.matrix = data[, c("GROUP.NUM", "")]
  #write.table(data, file = file.path(DATADIR, "image_list.txt"), row.names = F, quote = F, col.names = F)
  DATA = merge(DATA, covars.table, by = c("SUBJECT", "TP"), all.x = T)
  
  View(DATA)  
  results = list()
  results = vbanalysis(
    IMAGING_FILE,
    OUTPUT_DIR,
    DATA,
    MASK_FILE,
    MYTEST,
    remove_outliers = F, 
    to_gifti = to_gifti, 
    flip = flip
  )
  results$data = DATA
  save(results, file = file.path(OUTPUT_DIR, "results.RData"))
  return(results)
}

# check results: freeview -f /usr/local/freesurfer/subjects/fsaverage6/surf/lh.inflated:overlay=INTERCEPT_coef.func.gii

# results.average = doit(file.path(DATADIR, analysis_type), 
#                          image.list, 
#                          testaverage, 
#                          'tests/average',
#                          #MASK = 'mask_4mm.nii.gz',  
#                          MASK = 'mask.nii.gz',  
#                          IMAGES_NAME = 'image_list.txt',
#                          #IMAGING_NAME = 'images_4mm.nii.gz')
#                          IMAGING_NAME = 'images.nii.gz',
#                          to_gifti = mysurf)


if (analysis_type %in% c('volume', 'surfR', 'surfL')) 
{
  
  # tests for activation maps
  results.quadratic = doit(file.path(DATADIR, analysis_type), 
                           image.list, 
                           testquadraticrun, 
                           'tests/quadratic',
                           #MASK = 'mask_4mm.nii.gz',  
                           MASK = 'mask.nii.gz',  
                           IMAGES_NAME = 'image_list.txt',
                           #IMAGING_NAME = 'images_4mm.nii.gz')
                           IMAGING_NAME = 'images.nii.gz',
                           to_gifti = mysurf)

  results.comparison = doit(file.path(DATADIR, analysis_type), 
                            image.list, 
                            modelcomparisonrun, 
                            'tests/comparison', 
                            #MASK = 'mask_4mm.nii.gz',  
                            MASK = 'mask.nii.gz',  
                            IMAGES_NAME = 'image_list.txt',
                            #IMAGING_NAME = 'images_4mm.nii.gz')
                            IMAGING_NAME = 'images.nii.gz',
                            to_gifti = mysurf)

  
} else {
  
  # tests for marker of pattern variability
  results.linear_mv = doit(file.path(DATADIR, analysis_type), 
                           image.list, 
                           testlinear_mv, 
                           'tests/linear_mv',
                           MASK = 'mask.nii.gz',  
                           IMAGES_NAME = 'image_list.txt',
                           IMAGING_NAME = 'images.nii.gz',
                           to_gifti = mysurf)
  
}

stophere

results.linear = doit(file.path(DATADIR, analysis_type), 
                      image.list, 
                      testlinearrun, 
                      'tests/linear', 
                      #MASK = 'mask_4mm.nii.gz',  
                      MASK = 'mask.nii.gz',  
                      IMAGES_NAME = 'image_list.txt',
                      #IMAGING_NAME = 'images_4mm.nii.gz')
                      IMAGING_NAME = 'images.nii.gz')

results.asymptotic = doit(DATADIR, 
                         image.list, 
                         testasymptoticrun, 
                         'tests/asymptotic', 
                         #MASK = 'mask_4mm.nii.gz',  
                         MASK = 'mask.nii.gz',  
                         IMAGES_NAME = 'image_list.txt',
                         #IMAGING_NAME = 'images_4mm.nii.gz')
                         IMAGING_NAME = 'images.nii.gz')


stophere
# lme analyses
results.average = doit(DATADIR, 
                       image.list, 
                       testaverage, 
                       'tests/average', 
                       #MASK = 'mask_4mm.nii.gz',  
                       MASK = 'mask.nii.gz',  
                       IMAGES_NAME = 'image_list.txt',
                       #IMAGING_NAME = 'images_4mm.nii.gz')
                       IMAGING_NAME = 'images.nii.gz')
stophere

#IMAGING_NAME = 'images.nii.gz')

coords = c(75, 59, 65) # somatosensory
coords = c(54, 36, 50) # PCC
coords = c(27, 48, 59) # parietal

setwd(DATADIR)

images = readNIfTI('images.nii.gz')

radius = 5
get_ts = function(coords, radius, images){
        
system(paste('fslmaths mask -roi', 
             coords[1], '1', 
             coords[2], '1',
             coords[3], '1', '0 1 point_mask'))

system(paste('fslmaths point_mask -kernel sphere', radius, '-fmean -bin sphere_mask'))           
system('rm point_mask.nii.gz')           

mask =  readNIfTI('sphere_mask.nii.gz')

ts = NULL
for (i in seq(dim(images)[4])) {
  ts[i] <- images[, , , i][mask > 0]
}

return(ts)
}

y = get_ts(coords, radius, images)
results = results.quadratic
#y = results$imaging.mat[, which.min(results$pvalues[, "GROUP_x_TRAINING_p"])]
#y = results$imaging.mat[, which.min(results$pvalues[, "GROUP_p"])]
y = rnorm(nrow(X))
X = results$data
X$y = y


theme_lh = function () { 
  theme_bw(base_size=20)
}

#modelcomparisonrun(X, y)
#model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT) + (1|CONFIGURATION), data = X)
model = lmer(y ~ 1  + CONFIGURATION + GROUP*(TRAINING + TRAINING.Q) + 
               (1 | SUBJECT/TP) + (TRAINING + TRAINING.Q|SUBJECT), data = X)

plot(X$TRAINING, X$y, col = ifelse(X$GROUP == "Experimental", 'red', 'blue') )

se = function(x) sd(x, na.rm = T)/sqrt(length(x))

X.mean = X %>% group_by(GROUP, TP) %>% summarise(sem = se(y), y = mean(y))
myplot = ggplot(X.mean, aes(x = TP, group = GROUP, col = GROUP, y = y, ymin  = y - sem, ymax = y + sem)) + 
  geom_line() +  geom_errorbar() +  theme_lh()
print(myplot)

myplot = ggplot(X.mean, aes(x = TP, group = GROUP, col = GROUP, y = y)) + geom_line() 
print(myplot)

