#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:31:00 2020
@author: benjamin.garzon@ki.se
"""

# remove bad cases

# back to LSA and sync with pipeline
MINACC = 1.0
num_cores = 5
MINCORRECT = 3
NTYPES = 4

from joblib import Parallel, delayed
import pandas as pd
import os
import numpy as np
from glob import glob
from roi_analysis_funcs import process, gather_results
import seaborn as sns
import matplotlib.pyplot as plt
    
# try if gridsearch improves
# project the labels to volume
overwrite_extract = False 
overwrite_scores = True
output_data = False

#suffix = 'mask' # within mask using Glasser parc
#labels_file = '/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/motor_roi_parcels.txt'
#labels_file = '/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/mask_parcels.txt'
#labels_dir = '/data/lv0/MotorSkill/labels/fsaverage'

suffix = 'tess' # whole brain using tessellation
labels_file = '/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/tessellation162_parcels.txt'
labels_dir = '/data/lv0/MotorSkill/labels/fsaverage/tessellation162'
effects_name = 'effects.nii.gz'
if __name__ == "__main__":
    
    WD = '/data/lv0/MotorSkill/'
    subjects_dir = os.path.join(WD, 'fmriprep', 'freesurfer')

    analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis/'
    runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])

    skipped_file = '/data/lv0/MotorSkill/QC/skipped_%s.csv'%(suffix)
    
    paths = glob(os.path.join(analysis_dir, 'sub-lue1101', 'ses-*', 
                                  effects_name))
    
    # open label file
    labels = np.genfromtxt(labels_file, dtype = 'str').tolist() #[:9]
    skipped = []
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

        # reorder in time to sync with data!        
        allsequences = allsequences.sort_values(by=['run', 'onset'])

        sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
        sequences.loc[:, 'iscorrect'] = sequences.accuracy >= MINACC
        ncorrect = np.sum(sequences.iscorrect)
        valid_types = sequences.groupby(['run', 'seq_type']).iscorrect.sum()
        valid_runs = (valid_types >= MINCORRECT).groupby('run').sum() == NTYPES
        valid_runs = valid_runs.index[valid_runs]

        print("Subject %s session %d. Total correct: %d of %d"%(subject, session, ncorrect, 
                                         len(sequences.accuracy)))
        if len(valid_runs) < 2:
            skipped.append([subject, session, ncorrect, len(valid_runs)])
            continue
        effects_file = os.path.join(analysis_dir, subject, 
                                 'ses-%d'%session, effects_name)

        results_file = os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                                'surf', 'roi_scores_%s.csv'%(suffix))
        
        roi_data_file = os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                                'surf', 'roi_data_%s.csv'%(suffix)) \
            if output_data else None
        
        ribbon_img = os.path.join(subjects_dir, subject, 'mri', 'ribbon.mgz')
        ribbon_mask = os.path.join(analysis_dir, subject, 'label', 
                                   'ribbon_mask.nii.gz')
        if not os.path.exists(ribbon_mask) or overwrite_extract:
            command = 'mri_extract_label %s 3 42 %s;'%(ribbon_img, ribbon_mask) + \
                      'fslmaths %s -bin %s -odt char'%(ribbon_mask, ribbon_mask)
            os.system(command)
            
        if overwrite_scores or not os.path.exists(results_file):
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
                                                      roi_data_file,
                                                      do_prewhitening = True,
                                                      overwrite_extract = overwrite_extract, 
                                                      permutate = False)
                                                     for label in labels)
            
            mycols = ['subject', 'session', 'hemi', 'label', 'resnorm',
            'alpha_trained', 
            'alpha_untrained', 
            'alpha_trained_PCM', 
            'alpha_untrained_PCM', 
            'crossnobis_trained',
            'crossnobis_untrained',
            'svm_acc', 'svm_cm_acc', 
            'svm_acc_trained', 
            'svm_acc_untrained',
            'mean_signal_trained',
            'mean_signal_untrained',
            'roi_size']
            results = pd.DataFrame(processed_list, columns = mycols)
            results.to_csv(results_file, index = False, float_format = '%.3f')
    
    results_all = gather_results(analysis_dir, suffix)
    
    results_all.to_csv(os.path.join(analysis_dir, 
                                    'surf', 
                                    'roi_scores_%s.csv'%(suffix)),
                                    index = False,
                                    float_format = '%.3f')
    plt.figure()    
    ax = sns.boxplot(x = 'session', 
                y = 'svm_acc', 
                hue = 'hemi', 
                data = results_all)
    ax.set_ylim([0, 1])
    plt.figure()
    ax = sns.boxplot(x = 'session', 
            y = 'alpha_untrained', 
            hue = 'hemi', 
            data = results_all)
    plt.figure()    
    ax = sns.boxplot(x = 'session', 
            y = 'crossnobis_untrained', 
            hue = 'hemi', 
            data = results_all)    #ax.set_ylim([-1, 1])
    ax.set_ylim([-1, 2])
    results_all.groupby(['session', 'label']).max().svm_acc.max()
    
    skipped_df = pd.DataFrame(skipped, columns = ['Subject', 'Session', 
                                                  'ncorrect', 'valid_runs'])
    skipped_df.to_csv(skipped_file, index = False)