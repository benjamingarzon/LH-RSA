#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:31:00 2020
@author: benjamin.garzon@ki.se
"""

# remove bad cases

# prewhitening...
# back to LSA and sync with pipeline
minacc = 1.0
num_cores = 10
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd
import os
import numpy as np
import pandas as pd
from glob import glob
from roi_analysis_funcs import project_label, second_moment, \
    fit_second_moment, process, gather_results
    
# project the labels to volume

if __name__ == "__main__":
    
    WD = '/data/lv0/MotorSkill/'
    labels_file = '/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/mask_parcels.txt'
    subjects_dir = os.path.join(WD, 'fmriprep', 'freesurfer')

    analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis/'
    labels_dir = '/data/lv0/MotorSkill/labels/fsaverage'
    runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])
    
    paths = glob(os.path.join(analysis_dir, 'sub-*', 'ses-*', 
                                  'effects.nii.gz'))
    
    # open label file
    labels = np.genfromtxt(labels_file, dtype = 'str').tolist() #[:9]
    
    # iterate across paths
    for path in paths:

        pathels = os.path.normpath(path).split(os.sep)
        subject = pathels[-3]
        session = int(pathels[-2].split('-')[1])
        sequences_file = os.path.join(WD, 'responses', subject, 
                                  'ses-%d'%session, 'sequences.csv')
        run_paths = glob( os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                           'run*'))
        runs = np.array([int(r[-1]) for r in run_paths])
        # find runs
        allsequences = pd.read_csv(sequences_file, sep = ' ')
        sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
        iscorrect = sequences.accuracy >= minacc
        ncorrect = np.sum(iscorrect)
        
        print("Subject %s session %d. Total correct: %d of %d"%(subject, session, ncorrect, 
                                         len(sequences.accuracy)))
        
        effects_file = os.path.join(analysis_dir, subject, 
                                 'ses-%d'%session, 'effects.nii.gz')
        residuals_file = os.path.join(analysis_dir, subject, 
                                 'ses-%d'%session, 'res4d.nii.gz')

        processed_list = Parallel(n_jobs = num_cores, 
                                  require = "sharedmem", 
                                  verbose = 0)(delayed(process) 
                                                     (label, 
                                                      subject, 
                                                      session,
                                                      subjects_dir, 
                                                      analysis_dir, 
                                                      labels_dir, 
                                                      sequences, 
                                                      effects_file,
                                                      permutate = False)
                                                     for label in labels)

        results_file = os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                                'surf', 'roi_scores.csv')
        mycols = ['subject', 'session', 'hemi', 'label', 'alpha', 'resnorm', 
                  'svm_acc', 'roi_size']
        results = pd.DataFrame(processed_list, columns = mycols)
        results.to_csv(results_file, index = False, float_format = '%.3f')
    
    results_all = gather_results(analysis_dir)
    
    results_all.to_csv(os.path.join(analysis_dir, 
                                    'surf', 
                                    'roi_scores.csv'),
                                    index = False,
                                    float_format = '%.3f')
    
    import seaborn as sns
    ax = sns.boxplot(x = 'session', 
                y = 'svm_acc', 
                hue = 'hemi', 
                data = results_all)
    ax.set_ylim([0, 1])
    results_all.groupby(['session', 'label']).max().svm_acc.max()