
library(neurobase)
library(lme4)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(dplyr)
library(multcomp)

theme_lh = theme_classic(base_size = 13, base_family = "Arial")
myPalette = c("red", "blue")

sem = function(x) sd(x, na.rm = T)/sqrt(length(x))

create_vol_rois = function(DATADIR, TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR){
  tests = NULL
  rois = NULL
  myplots = NULL
  j = 1  
  
  MASK_FILE = file.path(DATADIR, MASK_NAME)
  TEST = file.path(DATADIR, TESTDIR, paste(TESTNAME, 'nii.gz', sep = '.'))
  ROI_FILE = file.path(DATADIR, TESTDIR, paste(TESTNAME, 'mask_sphere', paste0(radius*2, 'mm') , sep = '-'))
  command = paste("./make_spherical_roi.sh", file.path(DATADIR, TESTDIR), TEST, DISTANCE, radius, ROI_FILE)
  print(command)
  system(command)

    # plot data
  load(file.path(DATADIR, TESTDIR, 'results.rda'))
  mask <- fast_readnii(MASK_FILE)
  roimask <- readNIfTI(ROI_FILE)
  for (myroi in seq(dim(roimask)[4])){
      title = paste('ROI', myroi, sep = '-')
      print("###############################################################################")
      print(title)
      print("###############################################################################")
      roi = roimask[, , , myroi]
      roi_indices = which(roi[mask>0]>0)
      imaging.mat = results$imaging.mat[, roi_indices, drop = F]
      browser()
      X = results$data[-results$excluded, ]
      X$y = rowMeans(imaging.mat)    
      
      model = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.Q) + (1 |SUBJECT), data = X)
      print(summary(model))
      print(sum(!is.na(X$y)))

      X = X %>% filter(SYSTEM == "Classic")  %>% group_by(SUBJECT) %>% dplyr::mutate(TP.baseline = min(TP), n = sum(!is.na(y))) %>% filter(TP.baseline == 1)%>% 
        dplyr::mutate(y.dem = (y - y[TP.baseline])/mean(y, na.rm = T)*100)
      print(table(paste(X$GROUP, X$TP)))
      
      YMIN = min(X$y.dem, na.rm = T) 
      YMAX = max(X$y.dem, na.rm = T)
      YMIN = YMIN + .1*(YMAX-YMIN)
      YMAX = YMAX - .1*(YMAX-YMIN)
        
      X.sem = X %>% group_by(TP, GROUP) %>% dplyr::summarise(y.mean = mean(y.dem, na.rm = T), y.sem = sem(y.dem))
      myplot = ggplot(NULL) +  geom_line(data = X, aes(x = TP, group = SUBJECT, col = GROUP, y = y.dem), alpha = 0.2) + geom_point(alpha = 0.2) +
      scale_colour_manual(values = myPalette) + theme_lh + theme(legend.position = "bottom") + ggtitle(title) + geom_hline(yintercept = 0, size = 0.2)
        
      myplots[[j]] = myplot +
      geom_line(data = X.sem, aes(x = TP, group = GROUP, col = GROUP, y = y.mean)) + 
      geom_point(data = X.sem, aes(x = TP, group = GROUP, col = GROUP, y = y.mean)) + 
      geom_errorbar(data = X.sem, aes(x = TP, ymax = y.mean + y.sem, ymin = y.mean - y.sem, group = GROUP, col = GROUP)) + ylim(YMIN, YMAX)
        
      j = j + 1
    
  }
  
  FIG_MAP = paste(gsub("/", "-", TESTDIR), TESTNAME, 'map', sep = '-')
  FIG_ROI = paste(gsub("/", "-", TESTDIR), TESTNAME, 'roi', sep = '-')

  if (!file.exists(file.path(FIGS_DIR, paste0(FIG_MAP, '.png')))){
    # plot maps
    command = paste("./show_surface.sh", FIGS_DIR,
                    tests['lh'],
                    tests['rh'],
                    THR, 1,
                    FIG_MAP)
    print(command)
    system(command)
  }

  # plot roi
  if (!file.exists(file.path(FIGS_DIR, paste0(FIG_ROI, '.png')))){
    command = paste("./show_surface.sh", FIGS_DIR,
                    rois['lh'],
                    rois['rh'],
                    0.1, length(myplots),
                    FIG_ROI)
    print(command)
    system(command)
  }
  
  FIG_DATA = paste(gsub("/", "-", TESTDIR), TESTNAME, 'data', sep = '-')
  #browser()
  ggsave(filename = file.path(FIGS_DIR, paste0(FIG_DATA, '.png')), 
         plot = ggarrange(plotlist = myplots, nrow = ceiling(length(myplots)/NCOL), ncol = min(NCOL, length(myplots))),
         width = 21, 
         height = 13.2
  )
  
  return(list(myplots = myplots, ROI_FILE = ROI_FILE))
  
}

create_surf_rois = function(DATADIR,TESTDIR, TESTNAME, DISTANCE, radius, MASK_NAME, THR, PRECOMP_ROI = NULL, wDEPTH = F) {
  tests = NULL
  rois = NULL
  myplots = NULL
  j = 1  
  for (hemi in c("lh", "rh")){
    MASK_FILE = file.path(DATADIR, paste(hemi, MASK_NAME, sep = '.'))
    DEST = paste(TESTDIR, hemi, sep = '.')
    TEST = file.path(DATADIR, DEST, paste(TESTNAME, 'func.gii', sep = '.'))
    SURFACE_FILE = SURFACE_FILES[hemi]  
    if (is.null(PRECOMP_ROI)){
    ROI_FILE = file.path(DATADIR, DEST, paste(TESTNAME, 'mask_circle', paste0(radius*2, 'mm') , sep = '-'))
    command = paste("./make_circular_roi.sh", file.path(DATADIR, DEST), SURFACE_FILE, TEST, DISTANCE, radius, ROI_FILE, THR)
    print(command)
    system(command)
    } else {
      ROI_FILE = PRECOMP_ROI
      if(length(grep(hemi, ROI_FILE))==0) next
    }
    rois[hemi] = paste0(ROI_FILE, '.all.func.gii')
    tests[hemi] = TEST
    # plot data
    load(file.path(DATADIR, DEST, 'results.rda'))
    mask <- fast_readnii(MASK_FILE)
    roimask <- readNIfTI(ROI_FILE)
    if (sum(roimask) == 0) next
    
    #browser()
    for (myroi in seq(dim(roimask)[4])){
      title = paste('ROI', hemi, myroi, sep = '-')
      print("###############################################################################")
      print(title)
      print("###############################################################################")
      roi = roimask[, , , myroi]
      roi_indices = which(roi[mask>0]>0)
      imaging.mat = results$imaging.mat[, roi_indices]

      X = results$data[-results$excluded, ]
      X$y = rowMeans(imaging.mat)    
      #browser()
      if (!wDEPTH) model = lmer(y ~ 1 + SYSTEM + GROUP*(TRAINING + TRAINING.Q) + (1 |SUBJECT), data = X)
      else model = lmer(y ~ 1 + SYSTEM + GROUP*DEPTH*(TRAINING + TRAINING.Q)  + (1|SUBJECT), data = X)

            print(summary(model))
      print(sum(!is.na(X$y)))
      
      c.2        = c(0, 0, 0, 0, 0,  1,  0)
      c.4        = c(0, 0, 0, 0,  0, 0,  1)
      
      #cont.mat = rbind(c.2, c.4)
      #glh = glht(model, linfct = cont.mat, alternative='greater') 
      #sumglh = summary(glh, test = Ftest())
      #print(sumglh)
      # remove baseline and select
      
      #X = X %>% filter(SYSTEM=='Classic')  
      #browser()
      K = ifelse(wDEPTH, 20, 5)
      X = X %>% group_by(SUBJECT) %>% dplyr::mutate(TP.baseline = min(TP), n = sum(!is.na(y))) %>% filter(TP.baseline == 1)%>% 
        dplyr::mutate(y.dem = (y - y[TP.baseline])/mean(y, na.rm = T)*100)
      print(table(paste(X$GROUP, X$TP)))
      
      if (!wDEPTH){     
        YMIN = min(X$y.dem, na.rm = T) 
        YMAX = max(X$y.dem, na.rm = T)
        YMIN = YMIN + .1*(YMAX-YMIN)
        YMAX = YMAX - .1*(YMAX-YMIN)
        
      X.sem = X %>% group_by(TP, GROUP) %>% dplyr::summarise(y.mean = mean(y.dem, na.rm = T), y.sem = sem(y.dem))
      myplot = ggplot(NULL) +  geom_line(data = X, aes(x = TP, group = SUBJECT, col = GROUP, y = y.dem), alpha = 0.2) + geom_point(alpha = 0.2) +
        scale_colour_manual(values = myPalette) + theme_lh + theme(legend.position = "bottom") + ggtitle(title) + geom_hline(yintercept = 0, size = 0.2)
      
      myplots[[j]] = myplot +
        geom_line(data = X.sem, aes(x = TP, group = GROUP, col = GROUP, y = y.mean)) + 
        geom_point(data = X.sem, aes(x = TP, group = GROUP, col = GROUP, y = y.mean)) + 
        geom_errorbar(data = X.sem, aes(x = TP, ymax = y.mean + y.sem, ymin = y.mean - y.sem, group = GROUP, col = GROUP)) + ylim(YMIN, YMAX)
      
      } else {
        YMIN = min(X$y, na.rm = T) 
        YMAX = max(X$y, na.rm = T)
        YMIN = YMIN + .1*(YMAX-YMIN)
        YMAX = YMAX - .1*(YMAX-YMIN)
        
        X.sem = X %>% group_by(TP, GROUP, DEPTH) %>% dplyr::summarise(y.mean = mean(y, na.rm = T), y.sem = sem(y)) #.dem
        
#        myplot = ggplot(NULL) +  geom_line(data = X, aes(x = TP, group = SUBJECT, col = GROUP, y = y.dem), alpha = 0.2) + geom_point(alpha = 0.2) +
#          scale_colour_manual(values = myPalette) + theme_lh + theme(legend.position = "bottom") 
        
        myplots[[j]] = ggplot(NULL) + geom_line(data = X.sem, aes(x = TP, group = GROUP, col = GROUP, y = y.mean)) + 
           geom_point(data = X.sem, aes(x = TP, group = GROUP, col = GROUP, y = y.mean)) + 
           geom_errorbar(data = X.sem, aes(x = TP, ymax = y.mean + y.sem, ymin = y.mean - y.sem, group = GROUP, col = GROUP)) + 
           ylim(YMIN, YMAX) + theme_lh + theme(legend.position = "bottom")  + ggtitle(title) + geom_hline(yintercept = 0, size = 0.2) + 
           facet_grid(. ~ DEPTH) +
           scale_colour_manual(values = myPalette) 
        
      }
      
      j = j + 1
    }
  }
  
  if (is.null(PRECOMP_ROI)){
  FIG_MAP = paste(gsub("/", "-", TESTDIR), TESTNAME, 'map', sep = '-')
  FIG_ROI = paste(gsub("/", "-", TESTDIR), TESTNAME, 'roi', sep = '-')
  
  if (!file.exists(file.path(FIGS_DIR, paste0(FIG_MAP, '.png')))){

    # plot maps
    command = paste("./show_surface.sh", FIGS_DIR, 
                    tests['lh'], 
                    tests['rh'], 
                    THR, 
                    1,
                    FIG_MAP,
                    inflated.lh,
                    inflated.rh)
    print(command)
    system(command)
  }
  
  # plot roi
  if (!file.exists(file.path(FIGS_DIR, paste0(FIG_ROI, '.png')))){
    command = paste("./show_surface.sh", FIGS_DIR, 
                    rois['lh'], 
                    rois['rh'], 
                    0.1, 
                    length(myplots),
                    FIG_ROI, 
                    inflated.lh,
                    inflated.rh)
    print(command)
    system(command)
  }
  }
  FIG_DATA = paste(gsub("/", "-", TESTDIR), TESTNAME, 'data', sep = '-')
  #browser()
  ggsave(filename = file.path(FIGS_DIR, paste0(FIG_DATA, '.png')), 
         plot = ggarrange(plotlist = myplots, nrow = ceiling(length(myplots)/NCOL), ncol = min(NCOL, length(myplots))),
         width = 21, 
         height = 13.2
         )
  
  return(list(myplots = myplots, ROI_FILE = ROI_FILE))
  
}
