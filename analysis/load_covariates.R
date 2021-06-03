###############################################################
# Load covariate data
###############################################################
RESPONSES_DIR= '~/Data/LeftHand/Lund1/responses/'
NREPS = 10 # reps for inputation

library(dplyr)
library(missForest)
library(ggplot2)

motion = read.table(file.path('~/Data/LeftHand/Lund1/fmriprep/motion.csv'))
motion = motion %>% dplyr::rename(SUBJECT = subject, TP = session, RUN = run, GROUP = group) %>% 
  group_by(SUBJECT, TP) %>% summarise(FD = mean(FD), DVARS = mean(DVARS))
motion$SUBJECT = as.character(sapply(motion$SUBJECT, function (x) substring(x , 5, 12)))
                              
trials.table = read.table(file.path(RESPONSES_DIR, "complete_trials_table.csv"), header = T, sep = ";")
demo.data = read.table(file.path(RESPONSES_DIR, "SubjectData.csv"), header = T, sep = ",")
system.table = read.table(file.path(RESPONSES_DIR, "ClassicMTX8.csv"), header = T, sep = ";") %>% 
  dplyr::rename(SUBJECT = Subject, TP = Session, system = Classic.MTX8)
schedules = read.table('~/Data/LeftHand/Lund1/responses/all_lue_schedule_tables.csv', header = T)

aseg = read.table('/home/share/MotorSkill/freesurfer/asegstats.csv', header = T) %>% 
  dplyr::rename(SCAN = Measure.volume, TIV = EstimatedTotalIntraCranialVol) %>% 
  dplyr::mutate(SUBJECT = substring(SCAN, 5, 11), TP = substring(SCAN, 13, 13)) %>%
  dplyr::select(SUBJECT, TP, TIV) %>% group_by(SUBJECT) %>% dplyr::mutate(medianTIV = median(TIV), sdTIV = sd(TIV)) 
aseg = merge(aseg, system.table, by = c("SUBJECT", "TP"))
print(ggplot(aseg, aes(x = reorder(SUBJECT, medianTIV), y = TIV)) + geom_boxplot())

myvars = c('physical_activity', 
         'sleep', 
         'cigarettes', 
         'liquid', 
         'coffee', 
         'alcohol')

# aggregate
covars.table = trials.table %>% filter(mod(sess_num, 6) == 0)%>% group_by(username, sess_num) %>% 
  summarise_at(myvars, mean) %>% 
  mutate(SESS_HOME = sess_num, SUBJECT = username, TP = SESS_HOME %/% 6 + 1) 

covars.table = merge(covars.table, system.table, by = c("SUBJECT", "TP"), all.y = T)
covars.table = merge(covars.table, 
                     demo.data,
                     by.x = c("SUBJECT"), 
                     by.y = c("StudyID"), all.x = T) %>% dplyr::select(SUBJECT, TP, gender, reasoning,
                                                            physical_activity, sleep, cigarettes, 
                                                            liquid, coffee, alcohol, system) 


incomplete = !complete.cases(covars.table)
#covars.table[, myvars] = missForest(covars.table[, myvars])$ximp
#covars.table = missForest(covars.table)$ximp

mysubsets = NULL
for (i in 1:NREPS){
  print(i)
  mysubs = base::sample(unique(covars.table$SUBJECT), 45)
  mysubset = subset(covars.table, SUBJECT %in% mysubs| incomplete == T)
  mysubset$SUBJECT = droplevels(mysubset$SUBJECT)
  result = missForest(mysubset)
  print(result$OOBerror)
  mysubset = result$ximp
  mysubsets = rbind(mysubsets, mysubset %>% dplyr::mutate(REP = i))
}

covars.table = merge(mysubsets %>% group_by(SUBJECT, TP) %>% dplyr::summarise_if(is.numeric, mean, na.rm = T), 
                     unique(mysubsets[c("SUBJECT", "TP", "gender", "system")]),
                            by = c("SUBJECT", "TP"))  %>% dplyr::select(-REP)
                     
covars.table[covars.table < 0] = 0
View(covars.table[incomplete, ])

covars.table = covars.table %>% rename_all(.funs = toupper)

covars.table = merge(covars.table, schedules[c("SUBJECT", "CONFIGURATION", "SCHEDULE_GROUP")], by = "SUBJECT")

# which are the subjects with MTX8?
MTX8_subjects = (covars.table %>% group_by(SUBJECT) %>% summarise(NSYSTEMS = sum(SYSTEM == "MTX8")) %>% filter(NSYSTEMS > 0))$SUBJECT


# createlinear/quadratic/asymptotic regressors

training = seq(0, 6)*6
training.quadratic = -(training - max(training)/2)^2 # center around midpoint
min.quad = min(training.quadratic)
training.quadratic = training.quadratic - min.quad
training.asymptotic = c(training.quadratic[1:4], rep(training.quadratic[4], 3))  

training = scale(training, center = T, scale = F)/10
training.quadratic = scale(training.quadratic, center = T, scale = F)/100
training.cubic = scale(training.quadratic, center = T, scale = F)/1000
training.asymptotic = scale(training.asymptotic, center = T, scale = F)/100 
training.level = c(-1, NA, 1, 1, 1, 1, NA)
training.level = c(-1, 1, 1, NA, NA, NA, NA)

plot(training, training.quadratic, type = "b", xlim = c(-2, 2))
points(training, training, type = "b", col = "green")
points(training, training.asymptotic, type = "b", col = "red")

covars.table$GROUP = "Intervention"
covars.table$GROUP[grep("lue.2", covars.table$SUBJECT)] = "Control"
covars.table$GROUP.NUM = ifelse(covars.table$GROUP == "Intervention", 1, 0)  
covars.table = covars.table %>%  mutate(GROUP = as.factor(GROUP))

covars.table$TRAINING = training[covars.table$TP]
covars.table$TRAINING.Q = training.quadratic[covars.table$TP]
covars.table$TRAINING.C = training.cubic[covars.table$TP]
covars.table$TRAINING.A = training.asymptotic[covars.table$TP]

# interactions
covars.table$TRAINING_x_GROUP = covars.table$TRAINING*covars.table$GROUP.NUM 
covars.table$TRAINING.Q_x_GROUP = covars.table$TRAINING.Q*covars.table$GROUP.NUM 

# add motion parameters
covars.table = merge(covars.table, motion, by = c("SUBJECT", "TP"))

# get an estimate of trial by trial variability and check how it changes
NTRIALS = 10

ok_trials = subset(trials.table, accuracy == 1)
unpaced_trials = subset(trials.table, paced == 0 & trial_type != "missed")
ok_unpaced_trials = subset(unpaced_trials, accuracy == 1)
last_trials = subset(ok_unpaced_trials, trial > NTRIALS) %>% 
  group_by(username, sess_num, true_sequence, seq_train, group) %>% 
  dplyr::summarise(meanMT = mean(MT), sdMT = sd(MT))  %>% filter(mod(sess_num, 6) == 0) %>% 
  mutate(SESS_HOME = sess_num, SUBJECT = username, TP = SESS_HOME %/% 6 + 1) 

last_trials.mean = last_trials %>% group_by(SUBJECT, TP, seq_train, group) %>% summarise(sdMT.log = mean(log(sdMT))) %>% arrange(SUBJECT, TP, seq_train)
last_trials.diff = last_trials.mean %>% group_by(SUBJECT, TP, group) %>% mutate(sdMT.log.untrained_trained = c(diff(sdMT.log), NA)) %>% na.omit

# compute difference
last_trials.diff = merge(last_trials.diff, covars.table, by = c("SUBJECT", "TP"))

model.variability = lmer(sdMT.log.untrained_trained ~ 1 + GROUP*TRAINING + GROUP*TRAINING.Q + CONFIGURATION +  (1 + TRAINING + TRAINING.Q|SUBJECT),
                         data = subset(last_trials.diff)) #, GROUP == "Experimental")) #  (1|CONFIGURATION)

print(summary(model.variability))

# average variability
theme_lh <- function () { 
  theme_bw(base_size=18, base_family="Avenir")
}
# %>% group_by() %>% summarise(sdMT = mean(sdMT.log.untrained_trained))

# untrained variance > trained variance
#plot.sdMT.diff = ggplot(last_trials.diff, aes(x = TP, y = sdMT.log.untrained_trained, col = GROUP, group = GROUP)) + geom_smooth(lwd=1) + geom_point(size = 3) + theme_lh()
#print(plot.sdMT.diff)

