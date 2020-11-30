rm(list=ls())
# plots histograms of thickness and T1 values

library(oro.nifti)
library(ggplot2)
library(reshape2)
library(dplyr)
#SUBJECTS_DIR= '/home/benjamin.garzon/Data/LeftHand/Lund1/freesurfer/results'
WD = '/home/share/MotorSkill/freesurfer/results'
SUBJECTS_DIR = WD 
MASK_DIR = '/home/benjamin.garzon/Data/LeftHand/Lund1/labels/fsaverage'

find_outliers = function(X){
  print("Finding outliers")
  mycols = colnames(X)[(grep("\\bV", colnames(X)))]
  l = length(mycols)
  
  for (k in seq(l)){
    for (subject in unique(X$subject)){
      m = mean(X[X$subject == subject, mycols[k]])
      s = sd(X[X$subject == subject, mycols[k]])
      local_outlier = (abs(X[X$subject == subject, mycols[k]]- m) > 3*s)
      X[local_outlier, mycols[k]] = NA
    }
  }
  
  return(X)
}

tomat = function(myfile, mylist, hemi, MASK_FILE = NULL){
  print(myfile)
  imaging = readNIfTI(myfile)
  
  shape <- dim(imaging)
  n_scans <- shape[4]
  print(shape)
  
  print('Reshaping data')  
  
  if (is.null(MASK_FILE)) 
  {
    n_voxels <- prod(shape[1:3])
    imaging.mat <- matrix(0, n_scans, n_voxels) 
    print(paste('Number of voxels:', n_voxels))  
    for (i in seq(n_scans)) {
      perc = i/n_scans*100
      if (perc %% 20 == 0) print(paste("Read", perc, "%"))
      imaging.mat[i, ] <- imaging[, , , i]
    }
  }
  else
  {   
    mask <- readNIfTI(MASK_FILE)
    n_voxels <- sum(mask > 0)
    imaging.mat <- matrix(0, n_scans, n_voxels) 
    print(paste('Voxels in mask:', n_voxels))  
    for (i in seq(n_scans)) {
      perc = i/n_scans*100
      if (perc %% 20 == 0) print(paste("Read", perc, "%"))
      imaging.mat[i, ] <- imaging[, , , i][mask > 0]
    }
  }
  
  rm(imaging)
  file_list = read.table(mylist)
  colnames(file_list) = c("subject", "session")
  colnames(imaging.mat) = paste0("V", seq(ncol(imaging.mat)))
  #  return(list(imaging.mat = imaging.mat, 
  #              file_list = file_list))
  mydata = as.data.frame(cbind(file_list, imaging.mat))
  mydata$hemi = hemi
  mydata$group = "Experimental"
  mydata$group[grep("lue.2", mydata$subject)] = "Control"
  
  #  mydata = find_outliers(mydata)
  return(mydata)
}

loadT1 = function(hemi, invdepth)
{
  data = tomat(file.path(SUBJECTS_DIR, paste0(hemi, '.T1-', invdepth, '.nii.gz')), 
               file.path(SUBJECTS_DIR, paste0(hemi, '.T1-', invdepth, '.txt')), 
               hemi, MASK_FILE = ifelse(hemi == 'lh', lh.mask, rh.mask))
  data = melt(data, id.vars = c("subject", "session", "hemi", "group"), variable.name = "voxel")
  data$depth = 1 - as.numeric(invdepth) # invert
  data$value = as.numeric(data$value)
  data$session = as.integer(data$session)
  data$cluster = paste(data$subject, data$session, data$depth)
  return(data)
}

loadthickness = function(hemi)
{
  data = tomat(file.path(SUBJECTS_DIR, paste0(hemi, '.thickness.nii.gz')), 
               file.path(SUBJECTS_DIR, paste0(hemi, '.thickness.txt')), 
               hemi, MASK_FILE = ifelse(hemi == 'lh', lh.mask, rh.mask))
  data = melt(data, id.vars = c("subject", "session", "hemi", "group"), variable.name = "voxel")
  data$value = as.numeric(data$value)
  data$session = as.integer(data$session)
  data$cluster = paste(data$subject, data$session, data$depth)
  return(data)
}

lh.mask = file.path(MASK_DIR, 'lh.G_precentral.nii.gz')
rh.mask = file.path(MASK_DIR, 'rh.G_precentral.nii.gz')
lh.mask = file.path(MASK_DIR, 'lh.mask.nii.gz')
rh.mask = file.path(MASK_DIR, 'rh.mask.nii.gz')

thicknessdata.melt = rbind(
  loadthickness('lh'),
  loadthickness('rh')
)
if (T)
  T1data.melt = rbind(
    loadT1('lh', '0.20'), 
    loadT1('lh', '0.40'),
    loadT1('lh', '0.60'),
    loadT1('lh', '0.80'),
    loadT1('rh', '0.20'), 
    loadT1('rh', '0.40'),
    loadT1('rh', '0.60'),
    loadT1('rh', '0.80')
  )

# compute ICC
# s1 = 1
# s2 = 2
# mysubject = "lue2102"
# mydata = thicknessdata.melt
# function(mydata, s1, s2){
#   for (mysubject in unique(mydata$subject)){
#     x = mydata$value[ mydata$subject == mysubject & mydata$session == s1]
#     y = mydata$value[ mydata$subject == mysubject & mydata$session == s2]
#   }
# return(cor(x, y))
# }
#   
# stophere
#myplot.T1 = ggplot(subset(T1data.melt, subject == "lue1101"), aes(x = value, group = session)) + 
  
# thickness 
subjects = unique(thicknessdata.melt$subject)

se = function(x) sd(x)/sqrt(length(x))
thickness.mean = thicknessdata.melt %>% group_by(subject, group, session, hemi) %>% summarise(mean = mean(value), sd = se(value)) %>% arrange(hemi, mean)

myplot.thickness = ggplot(thicknessdata.melt %>% filter(subject %in% sample(subjects, 6)), aes(x = value, group = cluster, color = as.factor(hemi))) + 
  geom_density() + 
  facet_wrap( ~ subject + hemi, ncol = 6) + 
  theme_classic() + 
  xlab('Cortical thickness') + 
  ylab('Frequency') 
print(myplot.thickness)

#myplot.thickness = ggplot(thicknessdata.melt, aes(x = value, group = cluster)) + 
#  geom_density() + 
#  facet_grid(subject ~ hemi) + 
#  theme_classic() + 
#  xlab('Cortical thickness') + 
#  ylab('Frequency') 
#print(myplot.thickness)

myplot.thickness = ggplot(thickness.mean, aes(x = session, y = mean, group = subject, col = subject)) + 
  geom_point(size = 2) + 
  geom_line() + 
  facet_grid(hemi ~ . ) + 
  theme_classic() + 
  xlab('Session') +  
  ylab('Cortical thickness') 
print(myplot.thickness)

myplot.thickness = ggplot(thickness.mean, aes(x = reorder(subject, mean), y = mean, group = subject, fill = subject)) + 
  geom_boxplot() + 
  facet_grid(hemi ~ . ) + 
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab('Subject') + 
  ylab('Cortical thickness') 
print(myplot.thickness)


# T1
T1.mean = T1data.melt %>% group_by(subject, group, session, hemi, depth) %>% summarise(mean = mean(value), sd = se(value)) %>% arrange(hemi, depth, mean)

# do some plots

myplot.T1 = ggplot(T1data.melt, aes(x = value, color = as.factor(depth), group = cluster)) + 
  geom_density() + 
  facet_wrap( ~ subject + hemi) + #facet_grid(subject ~ hemi) + 
  theme_classic() + 
  xlab('T1 value') + 
  ylab('Frequency') 
print(myplot.T1)

myplot.T1 = ggplot(T1.mean, aes(x = session, y = mean, group = subject, col = subject)) + 
  geom_point(size = 2) + 
  geom_line() + 
  facet_grid(hemi ~ as.factor(depth) ) + 
  theme_classic() + 
  xlab('Session') + 
  ylab('T1 value') 
print(myplot.T1)

myplot.T1 = ggplot(T1.mean, aes(x = reorder(subject, mean), y = mean, group = subject, fill = subject)) + 
  geom_boxplot() + 
  facet_grid(hemi ~ as.factor(depth) ) + 
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab('Subject') + 
  ylab('T1 value') 
print(myplot.T1)


##############################################
stophere
myplot.thickness = ggplot(subset(thicknessdata.melt, subject == "lue1101"), aes(x = value, group = session)) + 
  geom_density() + 
  facet_grid(subject ~ hemi) + 
  theme_classic() + 
  xlab('Cortical thickness') + 
  ylab('Frequency') + theme(legend.title = "Depth")
#+ 
#  scale_colour_gradientn(colours = terrain.colors(7))
print(myplot.thickness)

rhdata.1 = find_outliers(rhdata[c("subject", "session","V29763")])

myplot = ggplot(rhdata, aes(x = session, group = subject, col = subject, y = V29763)) + ylim(0, 4) + geom_line() + facet_grid(. ~ subject) 
print(myplot)
myplot = ggplot(rhdata, aes(x = subject, group = subject, col = subject, y = V29763)) + geom_boxplot() + ylab ('Cortical thickness (mm)') 
print(myplot)
myplot = ggplot(rhdata, aes(x = session, y = V29763, group = group, col = group)) + geom_smooth()
print(myplot)
myplot = ggplot(rhdata.1, aes(x = session, y = V29763, group = subject, col = group)) + geom_line()
print(myplot)

plot(colMeans(mydata[-c(1, 2)]), type = 'l')
heatmap(as.matrix(mydata[-c(1, 2)][, seq(1, 150000, 1000)]), Colv = NA, Rowv = NA, col = terrain.colors(256), scale = "column")




