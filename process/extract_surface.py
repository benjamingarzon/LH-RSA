# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from __future__ import division
import os, sys
from mvpa2.suite import *
from mvpa2.clfs.svm import LinearCSVMC
from mvpa2.base.hdf5 import h5save, h5load
import pandas as pd
from mvpa2.measures import rsa
import pylab as pl
from scipy.spatial.distance import squareform
from nilearn import plotting
from nilearn.surface import load_surf_data
from sklearn.manifold import MDS
from mvpa2.clfs.ridge import RidgeReg
from sklearn.decomposition import PCA
from utils import get_RDM_metric

from utils import plot_mtx, aggregate_matrix_normalize, aggregate_matrix
    

if __debug__:
    from mvpa2.base import debug
    debug.active += ["SVS", "SLC"]

# Define surface and volume data paths:

if False:
    surfpath = sys.argv[1] #'/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/freesurfer/sub-101/surf'
    datapath = sys.argv[2] #'/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/fmriprep/sub-101/func/'
    sequences_fn = sys.argv[3] #'/home/benjamin.garzon/Data/LeftHand/Lundpilot1/responses/sequences.csv'
    hemi = sys.argv[4] # 'rh'
    radius = float(sys.argv[5]) # 'rh'
    
else:    
    datapath = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/analysis/sub-103/ses-3'
    sequences_fn = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/responses/sub-103/ses-3/sequences.csv'
    label_fn = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/freesurfer/sub-103/label/rh.cortex.label'
    hemi = 'rh'
    radius = 15.0
    labels = pd.read_csv(label_fn, sep = '\s+', skiprows = 2, header = None)
    nvols = 385
    TR = 1.2
    duration_offset = 4.9#5.0
    NRUNS = 5
    
epi_fn = os.path.join(datapath, 'data.nii.gz')
effects_fn = os.path.join(datapath, 'effects.nii.gz')
surfpath = os.path.join(datapath, 'surf')

#load engine 
qe = h5load(os.path.join(datapath, 'surf', '%s-%.1f-qe.gzipped.hdf5'%(hemi, radius)))
mask = qe.voxsel.get_mask()

sequences = pd.read_csv(sequences_fn, sep = ' ')

# arrange data and timing
targets = []
accuracy = []
events = []
trials = []
blocks = []
time_coords = np.arange(nvols)*TR
chunks = np.concatenate([x*np.ones(nvols) for x in np.arange(NRUNS)])
trial_time = []
for run in np.unique(sequences.run):
    sequences_run = sequences.loc[ sequences.run == run, :]

    myevents = []
    mytrials = []
    
    trial_duration = np.diff(sequences_run.fixation)[0]
    # filter out events
    for index, row in sequences_run.iterrows():
        myevents.append({'duration': row.fixation_duration, 
                         'onset': row.fixation, 
                         'event': 'fixation', 
                         'trial': row.trial})

        myevents.append({'duration': row.duration +  duration_offset, 
                         'onset': row.onset, 
                         'event': 'execution', 
                         'accuracy': row.accuracy,
                         'trial': row.trial})

        mytrials.append({'duration': trial_duration, 
                         'onset': row.fixation, 
                         'trial': row.trial,
                         'block': row.block,
                         'target': str(row.seq_type),
                         'accuracy': row.accuracy})
 
    accuracy.extend(events2sample_attr(mytrials, time_coords, noinfolabel=0,
                                  onset_shift=0.0, condition_attr='accuracy'))

    events.extend(events2sample_attr(myevents, time_coords, noinfolabel='rest',
                                  onset_shift=0.0, condition_attr='event'))

    targets.extend(events2sample_attr(mytrials, time_coords, noinfolabel='none', 
                                 onset_shift=0.0, condition_attr='target'))

#    trials.extend(events2sample_attr(myevents, time_coords, noinfolabel='rest',
#                                  onset_shift=0.0, condition_attr='trial'))

    onsets = np.array(events2sample_attr(mytrials, time_coords, noinfolabel=-1,
                                  onset_shift=0.0, condition_attr='onset'))

    trials.extend(events2sample_attr(mytrials, time_coords, noinfolabel=-1,
                                  onset_shift=0.0, condition_attr='trial'))

    blocks.extend(events2sample_attr(mytrials, time_coords, noinfolabel=-1,
                                  onset_shift=0.0, condition_attr='block'))
    
    trial_time.extend(time_coords - onsets)

trial_time = np.array(trial_time)

fds = fmri_dataset(samples = epi_fn,
                   targets = targets, 
                   chunks = chunks,
                   mask = mask)

fds.sa.time_coords = trial_time
fds.sa['trials'] = trials
fds.sa['blocks'] = blocks
fds.sa['events'] = events
fds.sa['accuracy'] = accuracy
fds.sa['trial_time'] = trial_time
poly_detrend(fds, polyord=1, chunks_attr='chunks')
#fds = fds[fds.sa.events == 'execution']


fds_effects = fmri_dataset(samples = effects_fn,
                   targets =  sequences.seq_type, 
                   chunks = sequences.run,
                   mask = mask)
fds_effects.sa['accuracy'] = sequences.accuracy


if hemi == 'lh':
    vertices = {'premotor': 129336,
            'visual': 405,
            'control': 85199,# temporal
            'somatosensory': 68135}
else:
    vertices = {'visual': 1522, 
                'premotor': 3368,
                'somatosensory': 6298,
                'premotor_ventral': 29195,
                'control': 24572}


def get_vertex_data(vertex, label, fds, fds_effects):
#    label = 'premotor'
#    vertex = vertices[label]
    print("#######")
    print(label)
    # get time courses
    roi_indices = qe.voxsel.volgeom.lin2ijk(qe.voxsel[vertex])
    fd_indices = np.concatenate([ np.where([ np.array_equal(x, y) for x in fds.fa.voxel_indices]) for y in roi_indices ]).ravel()
    
    # clean accuracy
    # extract data and remove outside of valid trials 
    fds_effects = fds_effects[:, fd_indices]

    fdsz = fds[:, fd_indices]
    zscore(fdsz, chunks_attr='chunks', param_est=('events', ['rest']),
            dtype='float32')    

    nifti_mask = qe.voxsel.get_nifti_image_mask([vertex])
    nifti_mask.to_filename(os.path.join(datapath, 'results', '%s-roi.nii.gz'%(label)))
    
    return((fdsz, fds_effects))

vertex_data = {}
for region in vertices.keys():
    vertex_data[region] = get_vertex_data(vertices[region], region, fds, fds_effects)
    
    
def plot_vertex(label, vertex_data):

    fdsz = vertex_data[label][0]
    fds_effects = vertex_data[label][1]
    
    print("#######")
    print(label)
        
    # clean accuracy
    # extract data and remove outside of valid trials 
    NCOMPS = 10
    pca = PCA(n_components = NCOMPS, whiten = False)
    nels = len(np.unique(fds_effects.targets))
    meanRDMs = []
    meanRDMs_pca = []
    
    for chunkind, chunk in enumerate(np.unique(fds_effects.chunks)):

        fds_red = fds_effects[fds_effects.chunks == chunk ]
    #use new class
        fds_pca = pca.fit_transform(fds_effects.samples)
        fds_pca_red = fds_pca[fds_effects.chunks == chunk ]
        pl.plot(pca.explained_variance_ratio_)  
        pl.xlim((0, NCOMPS))
    
        dist_matrix = rsa.pdist(fds_red[fds_red.sa.accuracy == 1, :], metric='correlation')
        dist_matrix_pca = rsa.pdist(fds_pca_red[fds_red.sa.accuracy == 1, :], metric='correlation')
        meanRDM, within, between = get_RDM_metric(dist_matrix, 'correlation', fds_red[fds_red.sa.accuracy == 1, :].targets)
        meanRDM_pca, within_pca, between_pca = get_RDM_metric(dist_matrix_pca, 'correlation', fds_red[fds_red.sa.accuracy == 1, :].targets)
        meanRDMs.append(meanRDM)
        meanRDMs_pca.append(meanRDM_pca)
           
    meanRDMs = 1 - np.tanh(np.nanmean(np.dstack(meanRDMs), axis = 2))
    meanRDMs_pca = 1 - np.tanh(np.nanmean(np.dstack(meanRDMs_pca), axis = 2))

    pl.figure
    pl.subplot(1, 2, 1)
    pl.imshow(meanRDMs)
    pl.subplot(1, 2, 2)
    pl.imshow(meanRDMs_pca)
    pl.savefig(os.path.join(datapath, 'results', 'matrix-%s.png'%(label)))       

    RDM = meanRDMs_pca # select one
    
    colors = {1: 'b', 2: 'g', 3: 'r', 4:'k'}
        #    markers = ['o', 'v', '+', 's']
    mymat = RDM.copy()
    # as barplot
    within_values = np.diag(RDM)
    np.fill_diagonal(mymat, 0)
    between_values = np.sum(mymat, axis = 0)/(nels - 1)

    fig = pl.figure(figsize=(15, 10), dpi=300)
    ax = fig.gca()
    ind = np.unique(fds_effects.targets)    
    p1 = ax.bar(ind - 0.2, within_values, color = 'green', width=0.4)
    p2 = ax.bar(ind + 0.2, between_values, color = 'red', width=0.4)
    ax.set_xticks(ind) 
    ax.set_ylim((0, .7))
    fig.legend((p1[0], p2[0]), ('Within', 'Between'))            
    fig.savefig(os.path.join(datapath, 'results', 'barplot%s.png'%(label)))   

    pl.figure(figsize=(15, 10), dpi=300)

    for t, target in enumerate(np.unique(fds_red.targets)):
            sel = np.logical_and(fds_red.targets == target, fds_red.sa.accuracy == 1)   
#            sel = np.logical_and(fds_red.chunks == 2, sel)   
            std = np.std(fds_pca_red[sel , :], axis = 0) 
            mean = np.mean(fds_pca_red[sel , :], axis = 0) 
            pl.plot(mean, color = colors[target])   
            pl.plot(fds_pca_red[sel , :].T, color = colors[target], linewidth = 0.2)   
            pl.errorbar(x = np.arange(len(mean)), y = mean, yerr = std, fmt = 'o')
    pl.savefig(os.path.join(datapath, 'results', 'PCA-%s.png'%(label)))   
    
# ADD ALL
#    mtgs = mean_group_sample(['targets', 'chunks'])
#    fds_mean = mtgs(fds[fds.sa.trials > 0, fd_indices])
    
    X = fdsz.samples

    pl.figure(figsize=(15, 10), dpi=300)
    colors = {'1': 'b', '2': 'g', '3': 'r', '4':'k'}
#    for chunk in np.unique(fds.chunks):
#        pl.subplot(1, NRUNS, chunk + 1)
#        for trial in range(1, 6):
    for target in ['1', '2', '3', '4']:
        sel = fds.targets == target #np.logical_and(fds.chunks == chunk, fds.targets == target) 
        if np.sum(sel): 
            df = pd.DataFrame({'t': np.round(fds.sa.trial_time[sel], 1), 'x': np.mean(X[sel, :], axis = 1)})
            xx = np.linspace(np.min(df.t), np.max(df.t), 20)
            yy = np.interp(xx, df.t, df.x)
#                df = df.groupby('t').mean()
#                pl.plot(df.index, df.x, color = colors[target])
            pl.plot(xx, yy, color = colors[target])
            pl.xlim((0, 25))
            pl.ylim((-1, 2))
            pl.axvline(x=1.0, linestyle = '--', color = 'k') # fixation   
            pl.axvline(x=3.1, linestyle = '--', color = 'k') # execution  
            pl.axvline(x=6.6, linestyle = '--', color = 'k') # end of execution 
            pl.axvline(x=6.6 + 0.5 + 6.7, linestyle = '--', color = 'r') # mean for next trial  
            pl.axhline(y=0.0, linestyle = '--', color = 'g')   
                    
    pl.savefig(os.path.join(datapath, 'results', 'timecourses-%s.png'%(label)))       

#########
    if True:
        svm = LinearCSVMC()    
    #    ridge = RidgeReg()
        
        NPERMS = 100
        permutator = AttributePermutator('targets', count=NPERMS)
        partitioner = NFoldPartitioner()
    
        distr_est = MCNullDist(permutator, tail='left', enable_ca=['dist_samples'])
        
        cv_svm = CrossValidation(svm, partitioner,
                             errorfx=lambda p, t: np.mean(p == t),
                             enable_ca=['stats'])
        
    
    #    cv_ridge = CrossValidation(ridge, partitioner,
    #                         errorfx=lambda p, t: np.mean(p == t),
    #                         enable_ca=['stats'])
    
        results_svm = cv_svm(fds_effects[fds_effects.sa.accuracy == 1, :]) 
    #    results_ridge = cv_ridge(fds) 
    
        cv_mc = CrossValidation(svm,
                             partitioner,
                             postproc=mean_sample(),
                             null_dist=distr_est,
                             enable_ca=['stats'])
    
#        cv_mc = CrossValidation(svm, partitioner,
#                         errorfx=mean_mismatch_error,
#                         null_dist=distr_est,
#                         enable_ca=['stats'])
        # run
        
        results_clf = cv_mc(fds_effects[fds_effects.sa.accuracy == 1, :]) 
        p = cv_mc.ca.null_prob
        print 'Accuracy:', np.mean(results_svm)
        print 'CV-errors:', np.ravel(results_clf)
        print 'Corresponding p-values:',  np.ravel(p)
    
#        pl.hist(fds.samples)

        
    if False:
        #compute matrix and aggregate across trials
        dist_matrix = squareform(rsa.pdist(fds_effects[fds_effects.sa.accuracy == 1, :], metric='correlation'))
 #       dist_matrix = squareform(rsa.pdist(fds, metric='mahalanobis'))
        meanRDM, elements, score, within, between = aggregate_matrix(dist_matrix, fds_effects.chunks, fds_effects.targets)    
        
        print(label, score, within, between)
#        pl.figure(figsize=(15, 9), dpi=300)
#        pl.plot(within, between)
#        pl.xlim((0, 2))
#        pl.ylim((0, 2))
#        pl.savefig(os.path.join(datapath, 'results', 'score-%s.png'%(label)))            
        
        pl.figure(figsize=(15, 9), dpi=300)
        
        for chunkind, chunk in enumerate(np.unique(fds_effects.chunks)):
            pl.subplot(1, 5, chunkind + 1)
            mtx = meanRDM[:, :, chunkind]
            pl.imshow(mtx, interpolation='nearest')
            pl.xticks(range(len(mtx)), elements, rotation=-45)
            pl.yticks(range(len(mtx)), elements)
            pl.title('Correlation distances')
            pl.clim((0, np.nanmax(mtx)))
#            if chunkind == 4:
#                pl.colorbar()
        
        pl.savefig(os.path.join(datapath, 'results', 'target_distances-%s.png'%(label)))            
    
        # MDS
        embedding = MDS(n_components=2, dissimilarity = 'precomputed')
        X_transformed = embedding.fit_transform(dist_matrix)
        colors = ['k', 'r', 'g', 'b']
        markers = ['o', 'v', 'P', 's']
        pl.figure(figsize=(27, 9), dpi=300)
        
        for c,  chunk in enumerate(np.unique(fds_effects.chunks)):
            pl.subplot(1, 5, c + 1)
            xmin = np.min(X_transformed[fds_effects.chunks == chunk, 0])*1.1
            xmax = np.max(X_transformed[fds_effects.chunks == chunk, 0])*1.1
            ymin = np.min(X_transformed[fds_effects.chunks == chunk, 1])*1.1
            ymax = np.max(X_transformed[fds_effects.chunks == chunk, 1])*1.1

            for t, target in enumerate(np.unique(fds_effects.targets)):
                sel = np.logical_and(fds_effects.chunks == chunk, fds_effects.targets == target)
                
                pl.scatter(X_transformed[sel , 0],  X_transformed[sel, 1], 
                           marker = markers[t], color = colors[t], s = 5, alpha = 0.8)   
                pl.xlim((xmin, xmax))
                pl.ylim((ymin, ymax))
                pl.scatter(np.mean(X_transformed[sel , 0], axis = 0),  
                           np.mean(X_transformed[sel, 1], axis = 0),  
                           marker = markers[t], color = colors[t], s = 50)   
        pl.savefig(os.path.join(datapath, 'results', 'MDS-%s.png'%(label)))            

        # show centers and radius

#        from sklearn.covariance import LedoitWolf

#        for chunk in [1, 2, 3, 4]:
#            print(chunk)
#            fds_chunk = fds[fds.chunks == chunk]
#            w = LedoitWolf().fit(fds_chunk.samples)
#            dist_matrix = squareform(rsa.pdist(fds_chunk, metric='mahalanobis', VI = w.covariance_))
#            embedding = MDS(n_components=2, dissimilarity = 'precomputed')
#            X_transformed = embedding.fit_transform(dist_matrix)
#            lim = np.max(np.abs(1.1*X_transformed))
#            colors = ['k', 'r', 'g', 'b']
#            markers = ['o', 'v', '+', 's']
#            pl.figure(figsize=(5, 3), dpi=100)
#            for t, target in enumerate(np.unique(fds_chunk.targets)):
#                sel = fds_chunk.targets == target
#                
#                pl.scatter(X_transformed[sel , 0],  X_transformed[sel, 1], marker = markers[t], color = colors[t])   
#                pl.xlim((-lim, lim))
#                pl.ylim((-lim, lim))
#

#vertices = {'premotor': 129336,
#            'visual': 38712,
#            'SMA': 134859,
#            'control': 85199,# temporal
#            'control2': 95936,# subcortical
#            'somatosensory': 68135}


    
#clf = LinearCSVMC()

#cv = CrossValidation(clf, NFoldPartitioner(),
#                     errorfx=lambda p, t: np.mean(p == t),
#                     enable_ca=['stats'])


for region in vertex_data.keys():
    plot_vertex(region, vertex_data)

stophere 
# Look at this pattern
#plot_mtx(squareform(sl_rsa_fds1.samples[:, index]),
#         fds1.sa.targets,
#         'Most consistent searchlight pattern correlation distances')

#plot_mtx(squareform(sl_cons_fds2.samples[:, index]),
#         fds2.sa.targets,
#         'Most consistent searchlight correlations', clim = (-1,1))

dist_path_fn = os.path.join(surfpath, '%s.sl_distinct%.1f.func.gii' % (hemi, radius))
cons_path_fn = os.path.join(surfpath, '%s.sl_consist%.1f.func.gii' % (hemi, radius))
spread_path_fn = os.path.join(surfpath, '%s.sl_spread_correlation_%.1f.func.gii' % (hemi, radius))


# plot the surface data
hemi_text= 'right' if hemi == 'rh' else 'left'
for view in ['medial', 'lateral']:
    plotting.plot_surf_stat_map(surf_mesh = os.path.join(surfpath, '%s.inflated' % (hemi)), 
                                stat_map = load_surf_data(dist_path_fn)[:, 1], 
                                bg_map = os.path.join(surfpath, '%s.sulc' % (hemi)), 
                                hemi = hemi_text,
                                title = 'Distinctiveness %s hemisphere'%(hemi_text),
                                view = view,
                                output_file = os.path.join(datapath, 'results', '%s.distinct.%s.png' % (hemi, view)),
                                alpha = 0.7, 
                                bg_on_data = True, 
                                colorbar = True, 
                                threshold = 0.2, 
                                cmap = 'cold_hot')
    plotting.plot_surf_stat_map(surf_mesh = os.path.join(surfpath, '%s.inflated' % (hemi)), 
                                stat_map = load_surf_data(cons_path_fn)[:, 1], 
                                bg_map = os.path.join(surfpath, '%s.sulc' % (hemi)), 
                                hemi = hemi_text,
                                title = 'Consistency %s hemisphere'%(hemi_text),
                                view = view,
                                output_file = os.path.join(datapath, 'results', '%s.consist.%s.png' % (hemi, view)),
                                alpha = 0.7, 
                                bg_on_data = True,
                                colorbar = True,
                                cmap = 'cold_hot')
    plotting.plot_surf_stat_map(surf_mesh = os.path.join(surfpath, '%s.inflated' % (hemi)), 
                                stat_map = load_surf_data(spread_path_fn)[:, 1], 
                                bg_map = os.path.join(surfpath, '%s.sulc' % (hemi)), 
                                hemi = hemi_text,
                                title = 'Inverse compactness %s hemisphere'%(hemi_text),
                                view = view,
                                output_file = os.path.join(datapath, 'results', '%s.spread_correlation.%s.png' % (hemi, view)),
                                alpha = 0.7, 
                                vmax = 1.3,
                                bg_on_data = True, 
                                colorbar = True, 
                                threshold = 1.1, 
                                cmap = 'cold_hot')

#from surfer import Brain
#brain = Brain("fsaverage", "lh", "inflated", subjects_dir = os.path.join(surfpath, '../..'))
#overlay_file = "example_data/lh.sig.nii.gz"
#brain.add_overlay(overlay_file)
#brain.overlays["sig"].remove()
#brain.add_overlay(overlay_file, min=5, max=20, sign="pos")

