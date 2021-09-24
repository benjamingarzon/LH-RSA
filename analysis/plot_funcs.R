
library(neurobase)
library(lme4)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(dplyr)
library(multcomp)
library(freesurfer)
library(stringr)

DPI = 500 # default

clean_measure_names = function(x, keep_first = F){
  
  if (keep_first)
  sapply(x, function(z) {
    z = gsub('acc_', '', z)
    str_to_title(paste(strsplit(as.character(z), '_')[[1]], collapse = ' '))})
  else
    sapply(x, function(z) {
      z = gsub('acc_', '', z)
      str_to_title(paste(strsplit(as.character(z), '_')[[1]][-1], collapse = ' '))})
  
  }

label_map = function(x){
  labels = list('SMA.label' = 'supplementary motor cortex',
                         'PS.label'  = 'primary sensorimotor cortex',
                         'PM.label'  = 'premotor cortex',
                         'SPL.label' = 'superior parietal lobule')
  if (x %in% names(labels)) return(labels[[x]]) else return(x)
} 

hemi_map = list('lh' = 'left',
                'rh' = 'right')

theme_lh = theme_classic(base_size = 13, base_family = "Arial") 
myPalette = c("red", "blue", "green", "yellow", "black")


markoutliersIQR = function(x){
  y = ifelse(x %in% boxplot(x, plot = F)$out, NA, x)
  return(y)
}
sem = function(x)
  sd(x, na.rm = T) / sqrt(length(x))

plot_activation_data = function(X, title, regressout = F, YMIN = 0, YMAX = 4) {
  X = X %>% mutate(WAVE = as.numeric(substring(SUBJECT, 4, 4)), y = y/100)# %>% filter (WAVE> 1) to percent signal change

  if(F){
  model = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + 
                 GROUP*CONDITION*(TRAINING + TRAINING.A) +
                 (1 + TRAINING + TRAINING.A|SUBJECT), data = X)

  print(summary(model.untrained))
  }
  
  model.untrained = lmer(y ~ 1 + FD + SYSTEM + CONFIGURATION + 
                 GROUP*(TRAINING + TRAINING.A) +
                 (1 + TRAINING + TRAINING.A|SUBJECT), data = subset(X, WAVE > 1 & CONDITION == 'UntrainedCorrect'))
  print(summary(model.untrained)$coefficients[c(7, 10, 11), ])
  
  print(sum(!is.na(X$y)))
  
  # regress out movement and system
  if (regressout) X$y = predict(model, X %>% mutate(FD = 0, 
                                                    SYSTEM = 'Classic',
                                                    CONFIGURATION = '1')) + resid(model) 
  
#  X = X %>% group_by(SUBJECT) %>% dplyr::mutate(TP.baseline = min(TP), n = sum(!is.na(y))) %>% filter(TP.baseline == 1) %>%
#    dplyr::mutate(y.dem = (y - y[TP.baseline]) / mean(y, na.rm = T) *                     100)
  print(table(paste(X$GROUP, X$TP)))
  
  #X = X %>% group_by(TP) %>% mutate(y = markoutliersIQR(y)) %>% filter(!is.na(y)) 

  #X.subject = X %>% group_by(SUBJECT, TP, GROUP, CONDITION) %>% 
  #  summarise( y = mean(y, na.rm = T)) %>% 
  #  mutate(SUBJCON = paste(SUBJECT, CONDITION)) %>% 
  #  group_by(SUBJECT, GROUP, CONDITION) %>% mutate(y.dem = y - mean(y)) %>% ungroup()
  
  X.sem = X %>% group_by(GROUP, TP, CONDITION) %>% 
    summarise(y.sem = sem(y), y.mean = mean(y, na.rm = T)) %>% 
    mutate(GROUP_CONDITION = paste(GROUP, CONDITION))
  
  #YMIN = min(X.sem$y.mean, na.rm = T)
  #YMAX = max(X.sem$y.mean, na.rm = T)
  #YMIN = YMIN - .3 * (YMAX - YMIN)
  #YMAX = YMAX + .3 * (YMAX - YMIN)
  # myplot = ggplot() +  geom_line(data = X.subject,
  #                                  aes(
  #                                    x = TP,
  #                                    group = SUBJCON,
  #                                    col = CONDITION,
  #                                    y = y.dem
  #                                  ),
  #                                  alpha = 0.2) + geom_point(alpha = 0.2) 
  
    X.sem$CONDITION = gsub('Correct', '', X.sem$CONDITION)
    myplot =  ggplot() + 
    geom_line(data = X.sem, aes(
      x = TP,
      group = GROUP_CONDITION, #GROUP
      col = CONDITION,
      linetype = GROUP,
      y = y.mean
    )) +
    geom_point(data = X.sem, aes(
      x = TP,
      group = GROUP_CONDITION,
      col = CONDITION,
      linetype = GROUP,
      y = y.mean
    )) +
    geom_errorbar(data = X.sem,
                  aes(
                    x = TP,
                    ymax = y.mean + y.sem,
                    ymin = y.mean - y.sem,
                    group = GROUP_CONDITION,
                    col = CONDITION,
                    linetype = GROUP
                  )) +
#    facet_grid(. ~ GROUP) +
    xlab('Test session') +
    ylab('% signal change') + 
    scale_colour_manual(values = myPalette) + theme_lh + 
    ggtitle(title) + ylim(YMIN, YMAX) +
    theme(legend.position = "bottom", 
          legend.title = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=18),
          axis.title=element_text(size=20, face="bold"))
    
    if (YMIN < 0 ) myplot = myplot + geom_hline(yintercept = 0, size = 0.2)  
    
  return(myplot)
}


plot_data = function(X, title, wDEPTH = F) {
  if (!wDEPTH) {
    model = lmer(y ~ 1 + SYSTEM + GROUP * (TRAINING + TRAINING.Q) + (1 |
                                                                       SUBJECT),
                 data = X)
    
    model.experimental = lmer(
      y ~ 1 + SYSTEM + (TRAINING + TRAINING.Q)  + (1 |
                                                     SUBJECT),
      data = subset(X, GROUP == "Experimental")
    )
    print("Only experimental: ")
    print(summary(model.experimental))
    
  }
  else {
    model = lmer(y ~ 1 + SYSTEM + GROUP * DEPTH * (TRAINING + TRAINING.Q)  + (1 |
                                                                                SUBJECT),
                 data = X)
  }
  print(summary(model))
  print(sum(!is.na(X$y)))
  
  #X = X %>% filter(SYSTEM=='Classic')
  #browser()
  K = ifelse(wDEPTH, 20, 5)
  X = X %>% group_by(SUBJECT) %>% dplyr::mutate(TP.baseline = min(TP), n = sum(!is.na(y))) %>% filter(TP.baseline == 1) %>%
    dplyr::mutate(y.dem = (y - y[TP.baseline]) / mean(y, na.rm = T) *
                    100)
  print(table(paste(X$GROUP, X$TP)))
  
  if (!wDEPTH) {
    YMIN = min(X$y.dem, na.rm = T)
    YMAX = max(X$y.dem, na.rm = T)
    YMIN = YMIN + .1 * (YMAX - YMIN)
    YMAX = YMAX - .1 * (YMAX - YMIN)
    
    X.sem = X %>% group_by(TP, GROUP) %>% dplyr::summarise(y.mean = mean(y.dem, na.rm = T),
                                                           y.sem = sem(y.dem))
    # myplot = ggplot(NULL) +  geom_line(data = X,
    #                                    aes(
    #                                      x = TP,
    #                                      group = SUBJECT,
    #                                      col = GROUP,
    #                                      y = y.dem
    #                                    ),
    #                                    alpha = 0.2) + geom_point(alpha = 0.2)
    myplot =  ggplot(NULL) + #myplot +
      geom_line(data = X.sem, aes(
        x = TP,
        group = GROUP,
        col = GROUP,
        y = y.mean
      )) +
      geom_point(data = X.sem, aes(
        x = TP,
        group = GROUP,
        col = GROUP,
        y = y.mean
      )) +
      geom_errorbar(data = X.sem,
                    aes(
                      x = TP,
                      ymax = y.mean + y.sem,
                      ymin = y.mean - y.sem,
                      group = GROUP,
                      col = GROUP
                    )) +
      scale_colour_manual(values = myPalette) + theme_lh + theme(legend.position = "bottom") + 
      ggtitle(title) + geom_hline(yintercept = 0, size = 0.2)
    #+ ylim(YMIN, YMAX)
    
  } else {
    YMIN = min(X$y, na.rm = T)
    YMAX = max(X$y, na.rm = T)
    YMIN = YMIN + .1 * (YMAX - YMIN)
    YMAX = YMAX - .1 * (YMAX - YMIN)
    
    X.sem = X %>% group_by(TP, GROUP, DEPTH) %>% dplyr::summarise(y.mean = mean(y, na.rm = T), y.sem = sem(y)) #.dem
    
    #        myplot = ggplot(NULL) +  geom_line(data = X, aes(x = TP, group = SUBJECT, col = GROUP, y = y.dem), alpha = 0.2) + geom_point(alpha = 0.2) +
    #          scale_colour_manual(values = myPalette) + theme_lh + theme(legend.position = "bottom")
    
    myplot = ggplot(NULL) + geom_line(data = X.sem, aes(
      x = TP,
      group = GROUP,
      col = GROUP,
      y = y.mean
    )) +
      geom_point(data = X.sem, aes(
        x = TP,
        group = GROUP,
        col = GROUP,
        y = y.mean
      )) +
      geom_errorbar(data = X.sem,
                    aes(
                      x = TP,
                      ymax = y.mean + y.sem,
                      ymin = y.mean - y.sem,
                      group = GROUP,
                      col = GROUP
                    )) +
      ylim(YMIN, YMAX) + theme_lh + theme(legend.position = "bottom")  + ggtitle(title) + geom_hline(yintercept = 0, size = 0.2) +
      facet_grid(. ~ DEPTH) +
      scale_colour_manual(values = myPalette)
    
  }
  return(myplot)
}

create_vol_rois = function(DATADIR,
                           TESTDIR,
                           TESTNAME,
                           DISTANCE,
                           radius,
                           MASK_NAME,
                           THR,
                           plot_function = plot_data) {
  tests = NULL
  rois = NULL
  myplots = NULL
  j = 1
  
  MASK_FILE = file.path(DATADIR, MASK_NAME)
  TEST = file.path(DATADIR, TESTDIR, paste(TESTNAME, 'nii.gz', sep = '.'))
  ROI_FILE = file.path(DATADIR,
                       TESTDIR,
                       paste(TESTNAME, 'mask_sphere', paste0(radius * 2, 'mm') , sep = '-'))

  command = paste(
    "./make_spherical_roi.sh",
    file.path(DATADIR, TESTDIR),
    TEST,
    DISTANCE,
    radius,
    ROI_FILE,
    THR
  )
  print(command)
  system(command)
  
  # plot data
  myfile = file.path(DATADIR, TESTDIR, 'results.rda')
  if (!exists(myfile)) myfile = file.path(DATADIR, TESTDIR, 'results.RData')
  load(myfile)
  mask <- fast_readnii(MASK_FILE)
  roimask <- readNIfTI(ROI_FILE)
  if (sum(roimask) == 0)
    return(NULL)
  
  for (myroi in seq(dim(roimask)[4])) {
    title = paste('ROI', myroi, sep = '-')
    print(
      "###############################################################################"
    )
    print(title)
    print(
      "###############################################################################"
    )
    
    roi = roimask[, , , myroi]
    roi_indices = which(roi[mask > 0] > 0)
    imaging.mat = results$imaging.mat[, roi_indices, drop = F]
    if (is.null(results$complete_data)) {
      X = results$data[-results$excluded,] }
    else {
      X = results$data
    }
    X$y = rowMeans(imaging.mat)
    myplots[[j]] = plot_function(X, title)
    j = j + 1
    
  }
  
  FIG_DATA = paste(gsub("/", "-", TESTDIR), TESTNAME, 'data', sep = '-')
  #browser()
  ggsave(
    filename = file.path(FIGS_DIR, paste0(FIG_DATA, '.png')),
    plot = ggarrange(
      plotlist = myplots,
      nrow = ceiling(length(myplots) / NCOL),
      ncol = min(NCOL, length(myplots))
    ),
    width = 21,
    height = 13.2
  )
  
  return(list(myplots = myplots, ROI_FILE = ROI_FILE))
  
}

create_surf_rois = function(DATADIR,
                            TESTDIR,
                            TESTNAME,
                            DISTANCE,
                            radius,
                            MASK_NAME,
                            THR,
                            MAXTHR, 
                            PRECOMP_ROI = NULL,
                            wDEPTH = F, 
                            plot_function = plot_data, annot = NULL, maxplots = 16) {
  tests = NULL
  rois = NULL
  myplots = NULL
  j = 1
  for (hemi in c("lh", "rh")) {
    MASK_FILE = file.path(DATADIR, paste(hemi, MASK_NAME, sep = '.'))
    DEST = paste(TESTDIR, hemi, sep = '.')
    TEST = file.path(DATADIR, DEST, paste(TESTNAME, 'func.gii', sep = '.'))
    SURFACE_FILE = SURFACE_FILES[hemi]
    if (is.null(PRECOMP_ROI)) {
      #ROI_FILE = file.path(DATADIR,
      #                     DEST,
      #                     paste(TESTNAME, 'mask_circle', paste0(radius * 2, 'mm') , sep = '-'))
      ROI_FILE = file.path(DATADIR,
                           DEST,
                           paste(TESTNAME, 'mask_cluster', sep = '-'))
      command = paste(
        "./make_cluster_rois.sh",
        file.path(DATADIR, DEST),
        SURFACE_FILE,
        TEST,
        DISTANCE,
        radius,
        ROI_FILE,
        THR,
        hemi
      )
      print(command)
      system(command)
    } else {
      ROI_FILE = PRECOMP_ROI
      if (length(grep(hemi, ROI_FILE)) == 0)
        next
    }
    rois[hemi] = paste0(ROI_FILE, '.all.func.gii')

    tests[hemi] = TEST
    # plot data 
    myfile = file.path(DATADIR, DEST, 'results.rda')
    if (!exists(myfile)) myfile = file.path(DATADIR, DEST, 'results.RData')
    
    if (file.exists(myfile)) {
      load(myfile)
    }
    else {
      print("No results file found.") 
      next
      }
    mask <- fast_readnii(MASK_FILE)
    roimask <- readNIfTI(ROI_FILE)
    if (sum(roimask) == 0)
      next
    
#    for (myroi in seq(dim(roimask)[4])) {
      for (myroi in seq(max(roimask))) {
        title = paste0('ROI ', hemi, '_', myroi)

#      roi = roimask[, , , myroi]
      roi = 1*(roimask == myroi)
      roi_indices = which(roi[mask > 0] > 0)
      
      if (!is.null(annot)) {
        annot_data = read_annotation(annot[[hemi]])
        masked_labels = annot_data$label[annot_data$label>0]
        mylabels = table(masked_labels[roi_indices])
        mylabels = mylabels/sum(mylabels)*100
        w = which(annot_data$colortable$code %in% names(mylabels[1]))
        region_name = paste(sapply(annot_data$colortable$label[w], label_map), collapse = '/')
        #title = paste0(title, ' (', hemi_map[[hemi]], ' ', region_name, ')')
        title = str_to_title(
          paste0(hemi_map[[hemi]], ' ', region_name)
        )
      }
      
      print(
        "###############################################################################"
      )
      print(title)
      print(
        "###############################################################################"
      )
      
      imaging.mat = results$imaging.mat[, roi_indices, drop = F]
      if (is.null(results$complete_data)) {
        X = results$data[-results$excluded,] }
      else {
        X = results$data
      }
      
      X$y = rowMeans(imaging.mat)
      myplots[[j]] = plot_function(X, title)
      j = j + 1
    }
  }

  if (length(myplots) == 0) return(NULL)
  #browser()
  if (is.null(PRECOMP_ROI)) {
    FIG_MAP = paste(gsub("/", "-", TESTDIR), TESTNAME, 'map', sep = '-')
    FIG_ROI = paste(gsub("/", "-", TESTDIR), TESTNAME, 'roi', sep = '-')
    
    if (!file.exists(file.path(FIGS_DIR, paste0(FIG_MAP, '.png')))) {
      if (!is.null(annot)) {
        annot.lh = paste0('annot=', annot[['lh']],':annot_outline=1')
        annot.rh = paste0('annot=', annot[['rh']],':annot_outline=1')
      } else
      { annot.lh = annot.rh = '' }
      # plot maps
      command = paste(
        "./show_surface.sh",
        FIGS_DIR,
        tests['lh'],
        tests['rh'],
        THR,
        MAXTHR,
        FIG_MAP,
        inflated.lh,
        inflated.rh,
        annot.lh, 
        annot.rh
      )
      print(command)
      system(command)
    }
    
    # plot roi
    if (!file.exists(file.path(FIGS_DIR, paste0(FIG_ROI, '.png')))) {
      command = paste(
        "./show_surface.sh",
        FIGS_DIR,
        rois['lh'],
        rois['rh'],
        0.1,
        length(myplots),
        FIG_ROI,
        inflated.lh,
        inflated.rh,
        annot.lh,
        annot.rh
      )
      print(command)
      system(command)
    }
  }
  FIG_DATA = paste(gsub("/", "-", TESTDIR), TESTNAME, 'data', sep = '-')
  NCOL = min(length(myplots), NCOLMAX)
  nplots = min(length(myplots), maxplots)
  print(NCOL)
  print(nplots)
    ggsave(
    filename = file.path(FIGS_DIR, paste0(FIG_DATA, '.png')),
    plot = ggarrange(
      plotlist = myplots[seq(nplots)],
      nrow = ceiling(nplots / NCOL),
      ncol = min(NCOL, length(myplots))
    ),
    width = NCOL*10,
    height = ceiling(nplots / NCOL)*8,
    limitsize = F,
    dpi = DPI
  )
  
  return(list(myplots = myplots, ROI_FILE = ROI_FILE))
  
}
