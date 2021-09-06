#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:31:00 2020
@author: benjamin.garzon@ki.se
"""

# permutation and no permutation
# prewhitening and no perm: session/run
# derivative /effect
# try with derivatives
#try w o prewhitening

MINACC = 1.0
MINCORRECT = 3
NTYPES = 4

from joblib import Parallel, delayed
import pandas as pd
import os
import argparse
import numpy as np
from glob import glob
from roi_analysis_funcs import process, gather_results
import seaborn as sns
import matplotlib.pyplot as plt
    

def do_analysis(WD, PERMUTATE, overwrite_extract, overwrite_scores, output_data, 
                do_prewhitening, suffix, labels_file, labels_dir, 
                effects_name, num_cores):
    
    subjects_dir = os.path.join(WD, 'fmriprep', 'freesurfer')

    analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis/'
    roi_data_dir = '/data/lv0/MotorSkill/fmriprep/analysis/roi_data/'
    runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])

    skipped_file = '/data/lv0/MotorSkill/QC/skipped_%s.csv'%(suffix)
    
    paths = glob(os.path.join(analysis_dir, 'sub-lue*', 'ses-*', 
                                  effects_name))
    # open label file
    labels = np.genfromtxt(labels_file, dtype = 'str').tolist()
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

        # reorder in time to sync  
        allsequences = allsequences.sort_values(by=['run', 'seq_type'])

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
        
        roi_data_file = os.path.join(roi_data_dir, 
                                     '%s_ses-%d_roi_data_%s'%(subject, session, suffix)) \
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
                                                      do_prewhitening = do_prewhitening,
                                                      overwrite_extract = overwrite_extract, 
                                                      permutate = PERMUTATE)
                                                     for label in labels)
            
            mycols = ['subject', 'session', 'hemi', 'label', 'resnorm',
            'alpha_trained', 
            'alpha_untrained', 
            'alpha_trained_PCM', 
            'alpha_untrained_PCM', 
            'G_hat_trained',
            'G_hat_untrained',
            'clf_acc', 
            'clf_acc_perm', 
            'clf_acc_trained', 
            'clf_acc_untrained',
            'clf_acc_trained_untrained',
            'xnobis_same', 
            'xnobis_different', 
            'xnobis_trained_same', 
            'xnobis_trained_different', 
            'xnobis_untrained_same', 
            'xnobis_untrained_different', 
            'xnobis_trained_untrained',
            'cosine_same', 
            'cosine_different', 
            'cosine_trained_same', 
            'cosine_trained_different', 
            'cosine_untrained_same', 
            'cosine_untrained_different', 
            'cosine_trained_untrained',
            'correlation_same', 
            'correlation_different', 
            'correlation_trained_same', 
            'correlation_trained_different', 
            'correlation_untrained_same', 
            'correlation_untrained_different', 
            'correlation_trained_untrained',
            'euclidean_same', 
            'euclidean_different', 
            'euclidean_trained_same', 
            'euclidean_trained_different', 
            'euclidean_untrained_same', 
            'euclidean_untrained_different', 
            'euclidean_trained_untrained',
            'mean_signal_trained',
            'mean_signal_untrained',
            'valid',
            'valid_runs',
            'ncolumns',
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
    ax = sns.boxplot(x = 'label', 
                y = 'clf_acc', 
                hue = 'hemi', 
                data = results_all)
    ax.set_ylim([0, 0.5])

    plt.figure()    
    ax = sns.boxplot(x = 'label', 
                y = 'clf_acc_perm', 
                hue = 'hemi', 
                data = results_all)
    ax.set_ylim([0, 0.5])

    plt.figure()    
    ax = sns.boxplot(x = 'session', 
            y = 'alpha_untrained', 
            hue = 'hemi', 
            data = results_all)
    
    plt.figure()    
    ax = sns.boxplot(x = 'session', 
            y = 'xnobis_trained_untrained', 
            hue = 'hemi', 
            data = results_all)    #ax.set_ylim([-1, 1])
    ax.set_ylim([-.1, .3])
    results_all.groupby(['session', 'label']).max().clf_acc.max()
    
    skipped_df = pd.DataFrame(skipped, columns = ['Subject', 'Session', 
                                                  'ncorrect', 'valid_runs'])
    skipped_df.to_csv(skipped_file, index = False)
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run surface-based multivariate roi analysis.')
    
    parser.add_argument('--permutate', dest='permutate', action='store_true')
    parser.add_argument('--overwrite_extract', dest='overwrite_extract', action='store_true')
    parser.add_argument('--overwrite_scores', dest='overwrite_scores', action='store_true')
    parser.add_argument('--output_data', dest='output_data', action='store_true')
    parser.add_argument('--do_prewhitening', dest='do_prewhitening', action='store')
    parser.add_argument('--suffix', dest='suffix', action='store')
    parser.add_argument('--labels_file', dest='labels_file', action='store')
    parser.add_argument('--labels_dir', dest='labels_dir', action='store')
    parser.add_argument('--effects_name', dest='effects_name', action='store')
    parser.add_argument('--WD', dest='WD', action='store')
    parser.add_argument('--num_cores', dest='num_cores', action='store', 
                        type=int)
    
    args = parser.parse_args()

    do_analysis(args.WD, args.permutate, args.overwrite_extract, 
                args.overwrite_scores, args.output_data, 
                args.do_prewhitening, args.suffix, args.labels_file, 
                args.labels_dir, args.effects_name, args.num_cores)