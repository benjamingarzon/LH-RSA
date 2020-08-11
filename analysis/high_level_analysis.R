rm(list = ls())
library(pracma)
library(ggplot2)
library(dplyr)

source('~/Software/ImageLMMR/ImageLMMR.R')
source('~/Software/LeftHand/analysis/mytests.R')
nifti_convert = '/home/share/Software/HCP/workbench/bin_rh_linux64/wb_command -metric-convert -to-nifti'

# take run into account
# define vars
WD = '~/Data/LeftHand/Lund1/fmriprep/analysis'

# mris_convert /usr/local/freesurfer/subjects/fsaverage6/surf/rh.white ~/Data/LeftHand/Lund1/fmriprep/analysis/surf/rh.fsaverage.white.ds.gii


schedules = read.table('~/Data/LeftHand/Lund1/responses/all_schedule_tables.csv', header = T)
contrastname = "UntrainedCorrect"
contrastnum = 2

contrastname = "UntrainedCorrect_TrainedCorrect"
contrastnum = 4

collect_data = T
analysis_type = 'surfL' # 'volume'
NPROCS = 20
# list files
contrast = ifelse(analysis_type == 'volume',
                  paste0("\\bcope", contrastnum, "\\.nii\\.gz$"),
                  paste0("\\bcope", contrastnum, "\\.func\\.gii$")
                  )

rh.gii = file.path(WD, 'surf', 'rh.fsaverage.white.ds.gii')
lh.gii = file.path(WD, 'surf', 'lh.fsaverage.white.ds.gii')

if (analysis_type == 'volume') mysurf = '' else mysurf = ifelse(analysis_type == "surfL", lh.gii, rh.gii)

image.list = list.files(WD, full.names = T, recursive = T, pattern = contrast)
image.list = image.list[grep(analysis_type, image.list)]

print(image.list)

training = seq(0, 6)*6
training.quadratic = -(training - max(training)/2)^2 
min.quad = min(training.quadratic)
training.quadratic = training.quadratic - min.quad
training.asymptotic = c(training.quadratic[1:4], rep(training.quadratic[4], 3))  

training = scale(training, center = T, scale = F)/10
training.quadratic = scale(training.quadratic, center = T, scale = F)/100
training.asymptotic = scale(training.asymptotic, center = T, scale = F)/100 
training.level = c(-1, NA, 1, 1, 1, 1, NA)
training.level = c(-1, 1, 1, NA, NA, NA, NA)

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
   setwd(WD)
   system("fslmerge -t all_masks `find .. | grep bold_space-MNI152NLin2009cAsym_brainmask`")
   system("fslmaths all_masks -Tmean -mul 100 higherlevel/mask_sum") #`fslnvols      higherlevel/mask_mean all_masks`
   system("fslmaths all_masks -Tmin higherlevel/mask_min")
   system("mri_vol2vol --mov /home/benjamin.garzon/Data/LeftHand/Lund1/cat12/mask.nii.gz --targ higherlevel/mask_min.nii.gz --out higherlevel/GMprob.nii.gz --regheader")
   system("fslmaths higherlevel/GMprob -thr 0.2 -bin -mul higherlevel/mask_min higherlevel/mask_min; mv mask_min higherlevel/volume/")
   system("fslmaths higherlevel/GMprob -thr 0.2 -bin -mul higherlevel/mask_sum higherlevel/mask_sum; mv mask_sum higherlevel/volume/")
   system("rm all_masks.nii.gz")
  
   setwd(DATADIR)  
   #system("rm mask.nii.gz")
   system("cd volume; ln -s mask_min.nii.gz mask.nii.gz")
   system("fslmerge -t images `cat image_list.txt`")
   system("flirt -in mask -ref mask -out mask_4mm -interp nearestneighbour -applyisoxfm 4")
   system("flirt -in images -ref mask_4mm -out images_4mm -interp nearestneighbour -applyisoxfm 4")
   
  } else {
    setwd(DATADIR)  
    #system("rm mask.nii.gz")
    
    for (filename in image.list) {
      print(filename)
      system(paste(nifti_convert, 
      filename, 
      paste0(filename, ".nii.gz")
      ))
    }
    
    system(sprintf("IM=`cat %s | sed 's/.gii/.gii.nii.gz/g'`; fslmerge -t %s/images $IM; rm $IM", 
                   image_list_file, analysis_type))
    
    system(sprintf("fslmaths %s/images -Tmean -bin %s/mask", 
                   analysis_type, analysis_type))
  }
}

extract_substr = function(x, pattern, start, stop){
  i = strfind(x, pattern)
  substr(x, i + start, i + stop)
} 

doit = function(WD, IMAGES, MYTEST, OD, 
                MASK = '', 
                IMAGES_NAME = '',
                IMAGING_NAME = '',
                to_gifti = ''){
  
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
  DATA = merge(DATA, schedules[c("SUBJECT", "CONFIGURATION", "SCHEDULE_GROUP")], by = "SUBJECT")
  DATA$GROUP = "Experimental"
  DATA$GROUP[grep("lue.2", DATA$SUBJECT)] = "Control"
  DATA$GROUP.NUM = ifelse(DATA$GROUP == "Experimental", 1, 0)  
  
  DATA$TRAINING = training[DATA$TP]
  DATA$TRAINING.Q = training.quadratic[DATA$TP]
  DATA$TRAINING.A = training.asymptotic[DATA$TP]

  #  DATA$TRAINING.L = training.level[DATA$TP]

    # interactions
  DATA$TRAINING_x_GROUP = DATA$TRAINING*DATA$GROUP.NUM 
  DATA$TRAINING.Q_x_GROUP = DATA$TRAINING.Q*DATA$GROUP.NUM 

  #design.matrix = data[, c("GROUP.NUM", "")]
  #write.table(data, file = file.path(DATADIR, "image_list.txt"), row.names = F, quote = F, col.names = F)
  
  View(DATA)  
  results = list()
  results = vbanalysis(
    IMAGING_FILE,
    OUTPUT_DIR,
    DATA,
    MASK_FILE,
    MYTEST,
    remove_outliers = F, to_gifti = to_gifti, flip = T
  )
  results$data = DATA
  save(results, file = file.path(OUTPUT_DIR, "results.RData"))
  return(results)
}


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
stophere
results.comparison = doit(DATADIR, 
                          image.list, 
                          modelcomparisonrun, 
                          'tests/comparison', 
                          #MASK = 'mask_4mm.nii.gz',  
                          MASK = 'mask.nii.gz',  
                          IMAGES_NAME = 'image_list.txt',
                          #IMAGING_NAME = 'images_4mm.nii.gz')
                          IMAGING_NAME = 'images.nii.gz')

results.linear = doit(DATADIR, 
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
results = results.linear
#y = results$imaging.mat[, which.min(results$pvalues[, "GROUP_x_TRAINING_p"])]
#y = results$imaging.mat[, which.min(results$pvalues[, "GROUP_p"])]
X = results$data
X$y = y


theme_lh = function () { 
  theme_bw(base_size=20)
}

#modelcomparisonrun(X, y)
model = lmer(y ~ 1 + GROUP*(TRAINING + TRAINING.Q) + (1 + TRAINING + TRAINING.Q|SUBJECT) + (1|CONFIGURATION), data = X)

plot(X$TRAINING, X$y, col = ifelse(X$GROUP == "Experimental", 'red', 'blue') )

se = function(x) sd(x, na.rm = T)/sqrt(length(x))

X.mean = X %>% group_by(GROUP, TP) %>% summarise(sem = se(y), y = mean(y))
myplot = ggplot(X.mean, aes(x = TP, group = GROUP, col = GROUP, y = y, ymin  = y - sem, ymax = y + sem)) + 
  geom_line() +  geom_errorbar() +  theme_lh()
print(myplot)

myplot = ggplot(X.mean, aes(x = TP, group = GROUP, col = GROUP, y = y)) + geom_line() 
print(myplot)

