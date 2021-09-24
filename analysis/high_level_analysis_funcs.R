

extract_substr = function(x, pattern, start, stop){
  i = strfind(x, pattern)
  if (!is.null(i)) substr = substr(x, i + start, i + stop)
  else substr = "-1"
  return(substr)
  
} 

check_empty_data = function(IMAGING_FILE){
  x = system(paste("fslstats -t", IMAGING_FILE, "-a -m"), intern = T)
  notempty = as.numeric(x) != 0
  return(notempty)
}


doit = function(WD, IMAGES, MYTEST, OD, 
                MASK = '', 
                IMAGES_NAME = '',
                IMAGING_NAME = '',
                conditions = NULL, motion = NULL, 
                to_gifti = '', flip = F, remove_outliers = F){
  setwd(WD)
  IMAGING_FILE = file.path(WD, IMAGING_NAME)
  MASK_FILE = file.path(WD, MASK)
  OUTPUT_DIR = file.path(WD, OD)
  
  # create regressors
  subject = sapply(IMAGES, extract_substr, "sub-", 4, 10)
  sess_num = sapply(IMAGES, extract_substr, "ses-", 4, 4) 
  run = sapply(IMAGES, extract_substr, "run", 3, 3)
  DATA = data.frame(NAME = unlist(IMAGES), SUBJECT = subject, 
                    TP = as.numeric(sess_num), 
                    RUN = as.numeric(run), CONDITION = conditions) 
  DATA$SUBJECT.NUM = as.numeric(as.factor(DATA$SUBJECT))
  #write.table(data, file = file.path(DATADIR, "image_list.txt"), row.names = F, quote = F, col.names = F)
  DATA = left_join(DATA, covars.table, by = c("SUBJECT", "TP"))
  DATA = left_join(DATA, motion, by = c("SUBJECT", "TP", "RUN"))  
  DATA$CONFIGURATION = as.factor(DATA$CONFIGURATION)
  DATA$WAVE = as.numeric(substring(DATA$SUBJECT, 4, 4))
  DATA$FD.clean = markoutliersIQR(DATA$FD)
  notempty = check_empty_data(IMAGING_FILE)
  DATA$notempty = notempty
  DATA.GROUPED = DATA %>% filter(!is.na(FD.clean) & notempty ) %>% 
    group_by(SUBJECT, TP) %>% 
    sample_n(1) %>% 
    group_by(SUBJECT) %>% 
    summarise(COUNT = n()) %>%filter(COUNT<4)
  
  excluded = which(DATA$SUBJECT %in% DATA.GROUPED$SUBJECT | is.na(DATA$FD.clean) | !DATA$notempty)
  print(paste("Excluding ", length(excluded), "observations, with these complete subjects:"))
  print(DATA.GROUPED$SUBJECT)
  View(DATA)  

  # replicate DATA for each condition
  results = list()
  results = vbanalysis(
    IMAGING_FILE,
    OUTPUT_DIR,
    DATA,
    MASK_FILE,
    MYTEST,
    remove_outliers = remove_outliers, 
    to_gifti = to_gifti, excluded = excluded, 
    flip = flip
  )
  results$complete_data = DATA
  save(results, file = file.path(OUTPUT_DIR, "results.RData"))
  print(paste("Saving results to ", file.path(OUTPUT_DIR, "results.RData")))
  for (f in list.files(OUTPUT_DIR, pattern = '*.gii', full.names = T)){
    system(paste('mri_convert', f, f))
  }
  return(results)
}
