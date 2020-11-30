#!/usr/bin/env python
import pandas as pd
import os, sys
import numpy as np

# December 2018
# benjamin.garzon@gmail.com 

# creates explanatory variables (EV) 
# For fMRI analysis
# EV1 = fixation
# EV2 = stretch
# EV3 = trained correct sequence
# EV4 = untrained correct sequence
# EV5 = trained incorrect sequence
# EV6 = untrained incorrect sequence

# discount dummy time

#INIT = 4 beats * 0.7 sec = 2.8 sec
#BUFFER_TIME = 0.5 sec
#DUMMY_TIME = 6 sec

FILENAME = sys.argv[1]
SUBJECT = sys.argv[2]
SESSION = int(sys.argv[3])
INIT_TIME = float(sys.argv[4])
BUFFER_TIME = float(sys.argv[5])
DUMMY_TIME = float(sys.argv[6])
STRETCH_DURATION = float(sys.argv[7])
SESSION_FACTOR = int(sys.argv[8])

OUTPUTNAME = 'EV'
AGG_OUTPUTNAME = 'AEV'
OUTPUTNAME_SINGLE = 'TRIAL'
OUTPUTNAME_OTHER = 'OTHER'

# Function def
def output_ev(EV, i, suffix, noext = False):
    filename = "run%d/%s"%(i, suffix) if noext else "run%d/%s.csv"%(i, suffix)
    
    EV.to_csv(filename, 
              header=False, columns=['onset', 'duration', 'value'], 
              index=False, 
              sep = ' ', 
              float_format='%.2f')

data = pd.read_csv(FILENAME, sep = ' ')

data = data.loc[data.username == SUBJECT, :]
data = data.loc[data.sess_num == (SESSION-1)*SESSION_FACTOR, :]

run = data['run'].values
runs = np.unique(run).astype(int)
labels, unique_labels = pd.factorize(data['true_sequence'])
labels = labels + 1
data['label'] = labels

# For representational analysis / decoding
sequences = pd.DataFrame({
    'subject': data['username'],
    'sess_num': data['sess_num'],
    'trial': data['trial'],
    'run' : data['run'],
    'seq_type' : labels,
    'seq_train': data['seq_train'],
    'onset': data['clock_execution'] + INIT_TIME - DUMMY_TIME, 
    'offset': data['clock_finished'] - BUFFER_TIME - DUMMY_TIME, 
    'duration': data['clock_finished'] - data['clock_execution'] - INIT_TIME - BUFFER_TIME, 
    'accuracy': data['accuracy'],
    'fixation': data['clock_fixation'] - DUMMY_TIME,
    'fixation_duration': data['clock_execution'] - data['clock_fixation'] + INIT_TIME,
    'block': data['block'],    
    'true_sequence': data['true_sequence'].str.replace(" ", "")

    }) 

sequences.to_csv("sequences.csv", header=True, index=False, sep = ' ', float_format='%.2f')
print("#Total correct: %d of %d"%(np.sum(sequences.accuracy == 1.0), len(sequences.accuracy)))
print("Type  Correct  Run")
for i in runs:
    # collect variables
    data_run = data.loc[data.run == i, :]
#    print("Run %d"%(i))
    counts = data_run.loc[data_run.accuracy == 1.].label.value_counts(sort=False)
    counts = counts.to_frame()
    counts['run'] = i 
    print(counts.to_string(header=False))
 
    fixation = data_run['clock_fixation'] - DUMMY_TIME
    fixation_duration = data_run['clock_execution'] - data_run['clock_fixation'] 
    wait = data_run['clock_execution'] - DUMMY_TIME
    onset = data_run['clock_execution'] + INIT_TIME - DUMMY_TIME 
    duration = data_run['clock_finished'] - data_run['clock_execution'] - INIT_TIME - BUFFER_TIME
    accuracy = data_run['accuracy']
    labels = data_run['label']
    seq_train = data_run['seq_train']
    stretch_trial = data_run['stretch']
    stretch_trial.values[-1] = 0 # last trial is irrelevant
    directory = 'run' + str(i)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # for fMRI
    execute = pd.DataFrame({'onset': onset, 'duration': duration, 'value': 1}) 
    correct = pd.DataFrame({'onset': onset[accuracy==1.0], 'duration': duration[accuracy==1.0], 'value': 1}) 
    incorrect = pd.DataFrame({'onset': onset[accuracy<1.0], 'duration': duration[accuracy<1.0], 'value': 1})
    
    trainedcorrect = pd.DataFrame({'onset': onset[(seq_train=='trained') & (accuracy==1.0)], 
                                   'duration': duration[(seq_train=='trained') & (accuracy==1.0) ], 
                                   'value': 1}) 
    untrainedcorrect = pd.DataFrame({'onset': onset[(seq_train=='untrained') & (accuracy==1.0)], 
                                     'duration': duration[(seq_train=='untrained') & (accuracy==1.0) ], 
                                     'value': 1}) 
    trainedincorrect = pd.DataFrame({'onset': onset[(seq_train=='trained') & (accuracy < 1.0) ], 
                                     'duration': duration[(seq_train=='trained') & (accuracy < 1.0)], 
                                     'value': 1}) 
    untrainedincorrect = pd.DataFrame({'onset': onset[(seq_train=='untrained') & (accuracy < 1.0)], 
                                       'duration': duration[(seq_train=='untrained') & (accuracy < 1.0)], 
                                       'value': 1}) 
    
    fixate = pd.DataFrame({'onset': fixation, 'duration': fixation_duration + INIT_TIME, 'value': 1}) 
    stretch = pd.DataFrame({'onset': fixation.iloc[ np.where(stretch_trial > 0)[0] + 1 ] - STRETCH_DURATION, 
                    'duration': STRETCH_DURATION, 
                    'value': 1}) 

# 3 regressors
    
#    output_ev(execute, i,  OUTPUTNAME + '1')
#    output_ev(fixate, i,  OUTPUTNAME + '2')
#    output_ev(stretch, i,  OUTPUTNAME + '3')

# 4 regressors
#    output_ev(correct, i, OUTPUTNAME + '1')
#    output_ev(incorrect, i, OUTPUTNAME + '2')
#    output_ev(fixate, i, OUTPUTNAME + '3')
#    output_ev(stretch, i, OUTPUTNAME + '4')

# 6 regressors
    output_ev(fixate, i, OUTPUTNAME + '1')
    output_ev(stretch, i, OUTPUTNAME + '2')    
    output_ev(trainedcorrect, i, OUTPUTNAME + '3')
    output_ev(untrainedcorrect, i, OUTPUTNAME + '4')
    output_ev(trainedincorrect, i, OUTPUTNAME + '5')
    output_ev(untrainedincorrect, i, OUTPUTNAME + '6')
    
    # Single trials for RSA analysis
    for j, xx in enumerate(onset):
        SINGLE = pd.DataFrame({'onset': onset.values[j], 'duration': duration.values[j], 'value': 1}, index = [0]) 
        SINGLE.to_csv("run%d/%s%d"%(i, OUTPUTNAME_SINGLE, j+1), header=False, columns=['onset', 'duration', 'value'], index=False, sep = ' ', float_format='%.2f')
        noj = np.arange(len(onset.values)) != j
        for label in np.unique(labels): 
            OTHER = pd.DataFrame({'onset': onset[np.logical_and(noj, label == labels) ], 'duration': duration[np.logical_and(noj, label == labels)], 'value': 1 })            
            OTHER.to_csv("run%d/%s%d_%d"%(i, OUTPUTNAME_OTHER, label, j+1), header=False, columns=['onset', 'duration', 'value'], index=False, sep = ' ', float_format='%.2f')
        
    # Add fixation and stretch
    output_ev(fixate, i, 'FIXATION', noext = True)
    output_ev(stretch, i, 'STRETCH', noext = True)

    # aggregate effects - for LSA
    output_ev(fixate, i, AGG_OUTPUTNAME + '1')
    output_ev(stretch, i, AGG_OUTPUTNAME + '2')    
    output_ev(incorrect, i, AGG_OUTPUTNAME + '3')    

    mylabels = labels.values
    for l in np.unique(mylabels):        
        correct = pd.DataFrame({'onset': onset[(accuracy==1.0) & (mylabels == l)], 
                                       'duration': duration[(accuracy==1.0) & (mylabels == l)], 
                                       'value': 1}) 

        output_ev(correct, i, AGG_OUTPUTNAME + 'cor%s'%(l))

