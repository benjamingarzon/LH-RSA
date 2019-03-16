#!/usr/bin/env python
import pandas as pd
import os, sys
import numpy as np

# December 2018
# benjamin.garzon@gmail.com 

# creates explanatory variables (EV) 
# For fMRI analysis
# EV1 = correct sequence
# EV2 = incorrect sequence
# EV3 = fixation
# EV4 = stretch

# discount dummy time

#INIT = 4 beats * 0.7 sec = 2.8 sec
#BUFFER_TIME = 0.5 sec
#DUMMY_TIME = 6 sec

FILENAME = sys.argv[1]
INIT_TIME = float(sys.argv[2])
BUFFER_TIME = float(sys.argv[3])
DUMMY_TIME = float(sys.argv[4])
STRETCH_DIFF = float(sys.argv[5])
STRETCH_DURATION = float(sys.argv[6])
OUTPUTNAME = 'EV'
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



data = pd.read_csv(FILENAME, sep = ';')

run = data['run'].values
runs = np.unique(run).astype(int)
labels, unique_labels = pd.factorize(data['true_sequence'])
data['label'] = labels
# For representational analysis / decoding
sequences = pd.DataFrame({
    'trial': data['trial'],
    'run' : data['run'],
    'seq_type' : labels + 1,
    'seq_train': data['seq_train'],
    'onset': data['clock_execution'] + INIT_TIME - DUMMY_TIME, 
    'offset': data['clock_finished'] - BUFFER_TIME - DUMMY_TIME, 
    'duration': data['clock_finished'] - data['clock_execution'] - INIT_TIME - BUFFER_TIME, 
    'accuracy': data['accuracy'],
    'fixation': data['clock_fixation'] - DUMMY_TIME,
    'fixation_duration': data['clock_execution'] - data['clock_fixation'] + INIT_TIME,
    'block': data['block']

    }) 

sequences.to_csv("sequences.csv", header=True, index=False, sep = ' ', float_format='%.2f')

for i in runs:
    # collect variables
    data_run = data.loc[data.run == i, :]
    fixation = data_run['clock_fixation'] - DUMMY_TIME
    fixation_duration = data_run['clock_execution'] - data_run['clock_fixation'] 
    wait = data_run['clock_execution'] - DUMMY_TIME
    onset = data_run['clock_execution'] + INIT_TIME - DUMMY_TIME 
    duration = data_run['clock_finished'] - data_run['clock_execution'] - INIT_TIME - BUFFER_TIME
    accuracy = data_run['accuracy']
    labels = data_run['label']
    
    directory = 'run' + str(i)
    if not os.path.exists(directory):
        os.makedirs(directory)

    # for fMRI
    execute = pd.DataFrame({'onset': onset, 'duration': duration, 'value': 1}) 
    correct = pd.DataFrame({'onset': onset[accuracy==1.0], 'duration': duration[accuracy==1.0], 'value': 1}) 
    incorrect = pd.DataFrame({'onset': onset[accuracy<1.0], 'duration': duration[accuracy<1.0], 'value': 1}) 
    fixate = pd.DataFrame({'onset': fixation, 'duration': fixation_duration + INIT_TIME, 'value': 1}) 
    stretch = pd.DataFrame({'onset': fixation.iloc[np.where(np.diff(fixation) > STRETCH_DIFF)[0] + 1] - STRETCH_DURATION, 
                    'duration': STRETCH_DURATION, 
                    'value': 1}) 

# 3 regressors
    
#    output_ev(execute, i,  OUTPUTNAME + '1')
#    output_ev(fixate, i,  OUTPUTNAME + '2')
#    output_ev(stretch, i,  OUTPUTNAME + '3')

# 4 regressors
    output_ev(correct, i,  OUTPUTNAME + '1')
    output_ev(incorrect, i,  OUTPUTNAME + '2')
    output_ev(fixate, i,  OUTPUTNAME + '3')
    output_ev(stretch, i,  OUTPUTNAME + '4')
    
    # Single trials for RSA analysis
    for j, xx in enumerate(onset):
        SINGLE = pd.DataFrame({'onset': onset.values[j], 'duration': duration.values[j], 'value': 1}, index = [0]) 
        SINGLE.to_csv("run%d/%s%d"%(i, OUTPUTNAME_SINGLE, j+1), header=False, columns=['onset', 'duration', 'value'], index=False, sep = ' ', float_format='%.2f')
        noj = np.arange(len(onset.values)) != j
        for label in unique_labels: 
            print(label)
            OTHER = pd.DataFrame({'onset': onset[np.logical_and(noj, label != labels) ], 'duration': duration[np.logical_and(noj, label !=  label)], 'value': 1 })            
            OTHER.to_csv("run%d/%s%d"%(i, OUTPUTNAME_OTHER, j+1), header=False, columns=['onset', 'duration', 'value'], index=False, sep = ' ', float_format='%.2f')
        
    # Add fixation and stretch
    output_ev(fixate, i,  'FIXATION', noext = True)
    output_ev(stretch, i,  'STRETCH', noext = True)
    