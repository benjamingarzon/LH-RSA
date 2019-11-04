#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 17:27:38 2019

@author: benjamin.garzon
"""
import numpy as np
import pandas as pd
import sys, os
from mvpa2.suite import *
from mvpa2.base.hdf5 import h5save, h5load
from utils import *
################################################################
# arguments
################################################################

if True:
    datapath = sys.argv[1] 
    sequences_fn = sys.argv[2] 
    subject = sys.argv[3] # 'rh'
    session = sys.argv[4] # 'rh'
    hemi = sys.argv[5] # 'rh'
    roi_fn = sys.argv[6]
    radius = float(sys.argv[7])
    NPROC = float(sys.argv[8])
    testdir = sys.argv[9]
    runs = np.array([int(x) for x in sys.argv[10].split(' ') ])
    pairwise_metric = sys.argv[11]
else:    
    datapath = '/home/benjamin.garzon/Data/LeftHand/Lund1/fmriprep/analysis/sub-lue1105/ses-1/'
    sequences_fn = '/home/benjamin.garzon/Data/LeftHand/Lund1/responses/sub-lue1105/ses-1/sequences.csv'
    subject = 'sub-lue1105'
    session = 1
    roi_fn = '/home/benjamin.garzon/Data/LeftHand/Lund1/fmriprep/analysis/sub-lue1105/surf/spread_correlation_maxima.csv'
    hemi = "rh"
    radius = 15.0
    NPROC = 10
    testdir = 'metrics'
    runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])
    pairwise_metric = 'correlation'
    
NCOMPS = 10 
filter_accuracy = True   
minacc = 1.0 
################################################################
# function definitions
################################################################

def get_vertex_data(vertex, label, fds_effects):

    print("#######")
    print(label)
    # get time courses
    roi_indices = qe.voxsel.volgeom.lin2ijk(qe.voxsel[vertex])
    fd_indices = np.concatenate([ np.where([ np.array_equal(x, y) for x in fds.fa.voxel_indices]) for y in roi_indices ]).ravel()
    
    fds_effects = fds_effects[:, fd_indices]

    nifti_mask = qe.voxsel.get_nifti_image_mask([vertex])
    nifti_mask.to_filename(os.path.join(datapath, 'surf', 'rois', '%s.%s-roi.nii.gz'%(hemi, label)))
    
    return(fds_effects)


################################################################


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

#roi_ids = np.intersect1d(labels[0].ravel(), qe.ids)

################################################################
# function definitions
################################################################
correct = sequences.loc[ sequences.accuracy >= minacc, ["seq_type", "seq_train", "true_sequence", "run"] ]
correct = correct.groupby(["seq_type", "seq_train", "true_sequence", "run"]).size().to_frame(name = 'count').reset_index()
correct = correct[correct["count"] > 1 ]

#correct[["seq_type", "seq_train", "true_sequence"]].drop_duplicates().sort_values("seq_type").to_csv(
#    os.path.join(surfpath, 'sequence_types.csv'), 
#    header= False, 
#    index= False, sep = ' ', 
#    float_format='%.2f')

#correct.sort_values("seq_type").to_csv(
#    os.path.join(surfpath, 'sequence_counts.csv'), 
#    header= False, 
#    index= False, sep = ' ', 
#    float_format='%.2f')

rois = rois.loc[ rois.HEMI == hemi ]
value_list = []            

for i, row in rois.iterrows():
    print(row)

    vertex_data = get_vertex_data(row.INDEXMAX, row.NAME, fds)
    
#        dsm_ic = Pwithin_spread(mapper = flatten_mapper(), 
#            NCOMPS = 10, 
#            filter_accuracy = True,
#            accuracy = fds.sa.accuracy,
#            pairwise_metric = 'correlation',
#            square=False)
    
#        results_rsa = dsm_ic(vertex_data)                        
    
    ds = vertex_data
    mapped_ds = flatten_mapper()(ds)
    data = mapped_ds.samples
    
    # remove constant columns
    constantcols = np.all(data[1:] == data[:-1], axis=0)
    data = data[:, ~constantcols]
    if np.sum(constantcols)>0:
        print(np.sum(constantcols))
    try:               
        pca = PCA(n_components = NCOMPS, whiten = False)
        data_pca = pca.fit_transform(data)
        
        if filter_accuracy:
            data_pca = data_pca[fds.sa.accuracy >= minacc]
            chunks = ds.chunks[fds.sa.accuracy >= minacc]
            targets = ds.targets[fds.sa.accuracy >= minacc ]
        else:
            chunks = ds.chunks
            targets = ds.targets
    #                print(data_pca.shape)
    
        unique_chunks = np.unique(chunks) 
        Nchunks = len(unique_chunks)

        for i, seq_info in enumerate(correct.iterrows()):
            chunk = seq_info[1].run

            ds_chunk = data_pca[chunks == chunk]
            dsm = pdist(ds_chunk, 
                    metric=pairwise_metric)
            within, between = get_spread_chunk(dsm, 
                               pairwise_metric, 
                               targets[chunks == chunk])
                                
            for j, value in enumerate(within[seq_info[1].seq_type][0]):                     
                value_list.append(
                        (subject, session, chunk,
                         row.HEMI, row.NAME, row.INDEXMAX,
                         seq_info[1].seq_type, 
                         seq_info[1].seq_train, 
                         seq_info[1].true_sequence, 
                         j + 1, 
                         value, 'within')
                        ) 

            for j, value in enumerate(between[seq_info[1].seq_type][0]):                     
                value_list.append(
                        (subject, session, chunk,
                         row.HEMI, row.NAME, row.INDEXMAX,
                         seq_info[1].seq_type, 
                         seq_info[1].seq_train, 
                         seq_info[1].true_sequence, 
                         j + 1, 
                         value, 'between')
                        ) 
    except LinAlgError:
        print("PCA did not converge!")
    except ValueError:
        print("Value error")

#columns
#to data path
results = pd.DataFrame(value_list)
#results.columns = ['subject', 'session', 'run', 'hemi', 'name', 'vertex_index', 'seq_type', 
#               'seq_train', 'true_sequence', 'number', 'value', 
#               'dist_type']
results.to_csv(
    os.path.join(surfpath, '%s.roi_%s_distances.csv'%(hemi, pairwise_metric)), 
    index= False, sep = ' ', header = False, 
    float_format='%.2f'
    )