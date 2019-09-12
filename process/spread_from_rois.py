#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 17:27:38 2019

@author: benjamin.garzon
"""
import numpy as np
import pandas as pd
import sys

if False:
    datapath = sys.argv[1] 
    sequences_fn = sys.argv[2] 
    label_fn = sys.argv[3]
    hemi = sys.argv[4] # 'rh'
    radius = float(sys.argv[5])
    NPROC = float(sys.argv[6])
    testdir = sys.argv[7]
    runs = np.array([int(x) for x in sys.argv[8].split(' ') ])
else:    
    datapath = '/home/benjamin.garzon/Data/LeftHand/Lund1/fmriprep/analysis/sub-lue1105/ses-1/'
    sequences_fn = '/home/benjamin.garzon/Data/LeftHand/Lund1/responses/sub-lue1105/ses-1/sequences.csv'
    roi_fn = '/home/benjamin.garzon/Data/LeftHand/Lund1/fmriprep/analysis/sub-lue1105/surf/spread_correlation_maxima.csv'
    hemi = 'rh'
    radius = 15.0
    NPROC = 10
    testdir = 'metrics'
    runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])

rois = pd.read_csv(roi_fn, sep = '\s+')

surfpath = os.path.join(datapath, 'surf', testdir)
epi_fn1 = os.path.join(datapath, 'effects.nii.gz')
epi_fn2 = os.path.join(datapath, 'derivatives.nii.gz')

#load engine 
qe = h5load(os.path.join(datapath, 'surf', '%s-%.1f-qe.gzipped.hdf5'%(hemi, radius)))
mask = qe.voxsel.get_mask()

allsequences = pd.read_csv(sequences_fn, sep = ' ')
sequences = allsequences.loc[ allsequences.run.isin(runs), : ]

print("Total correct: %d of %d"%(np.sum(sequences.accuracy >= minacc), len(sequences.accuracy)))


fds1 = fmri_dataset(samples = epi_fn1,
                   targets =  sequences.seq_type, 
                   chunks = sequences.run,
                   mask = mask)

fds2 = fmri_dataset(samples = epi_fn2,
                   targets = sequences.seq_type, 
                   chunks = sequences.run,
                   mask = mask)

fds1.samples = stats.zscore(fds1.samples, axis = None)
fds2.samples = stats.zscore(fds2.samples, axis = None)

fds = fds1.copy()
fds.samples = np.dstack((fds1.samples, fds2.samples))
fds.sa.accuracy = sequences.accuracy

roi_ids = np.intersect1d(labels[0].ravel(), qe.ids)


for i, roi in rois.iterrows():
    print(roi.NAME)
    