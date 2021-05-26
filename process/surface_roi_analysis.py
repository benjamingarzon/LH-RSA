#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:31:00 2020
@author: benjamin.garzon@ki.se
"""

# remove bad cases

# back to LSA and sync with pipeline
MINACC = 1.0
num_cores = 40
MINCORRECT = 3
NTYPES = 4

import multiprocessing
from joblib import Parallel, delayed
import pandas as pd
import os
import numpy as np
import pandas as pd
from glob import glob
from roi_analysis_funcs import project_label, second_moment, \
    fit_second_moment, process, gather_results
import seaborn as sns
    
# project the labels to volume
overwrite = False
if __name__ == "__main__":
    
    WD = '/data/lv0/MotorSkill/'
    labels_file = '/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/motor_roi_parcels.txt'
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
        # find runs to sync with effects
        allsequences = pd.read_csv(sequences_file, sep = ' ')
        sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
        sequences.loc[:, 'iscorrect'] = sequences.accuracy >= MINACC
        ncorrect = np.sum(sequences.iscorrect)
        valid_types = sequences.groupby(['run', 'seq_type']).iscorrect.sum()
        valid_runs = (valid_types >= MINCORRECT).groupby('run').sum() == NTYPES
        valid_runs = valid_runs.index[valid_runs]

        print("Subject %s session %d. Total correct: %d of %d"%(subject, session, ncorrect, 
                                         len(sequences.accuracy)))
        if len(valid_runs) < 2:
            continue
#        if subject == 'sub-lue5103': 
#            continue
        effects_file = os.path.join(analysis_dir, subject, 
                                 'ses-%d'%session, 'effects.nii.gz')

        results_file = os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                                'surf', 'roi_scores.csv')
        
        ribbon_img = os.path.join(subjects_dir, subject, 'mri', 'ribbon.mgz')
        ribbon_mask = os.path.join(analysis_dir, subject, 'label', 
                                   'ribbon_mask.nii.gz')
        if not os.path.exists(ribbon_mask) or overwrite:
            command = 'mri_extract_label %s 3 42 %s;'%(ribbon_img, ribbon_mask) + \
                      'fslmaths %s -bin %s -odt char'%(ribbon_mask, ribbon_mask)
            os.system(command)
                                                                      
        if overwrite or not os.path.exists(results_file):
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
                                                      do_prewhitening = True,
                                                      permutate = False)
                                                     for label in labels)

            mycols = ['subject', 'session', 'hemi', 'label', 'resnorm',
            'alpha_trained', 
            'alpha_untrained', 
            'alpha_trained_PCM', 
            'alpha_untrained_PCM', 
            'crossnobis_trained',
            'crossnobis_untrained',
            'svm_acc', 'roi_size']
            results = pd.DataFrame(processed_list, columns = mycols)
            results.to_csv(results_file, index = False, float_format = '%.3f')
    
    results_all = gather_results(analysis_dir)
    
    results_all.to_csv(os.path.join(analysis_dir, 
                                    'surf', 
                                    'roi_scores.csv'),
                                    index = False,
                                    float_format = '%.3f')
    
    ax = sns.boxplot(x = 'session', 
                y = 'svm_acc', 
                hue = 'hemi', 
                data = results_all)
    ax.set_ylim([0, 1])
    ax = sns.boxplot(x = 'session', 
            y = 'alpha_trained', 
            hue = 'hemi', 
            data = results_all)
    ax.set_ylim([-1, 1])

    results_all.groupby(['session', 'label']).max().svm_acc.max()