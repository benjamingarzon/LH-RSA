#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:34:38 2020

@author: benjamin.garzon@gmail.com
"""
# compute LSS from LSA effects

# Usage:compute_effects.py X_LSA_file beta_LSA_file res_LSA_file 
#sequences_fn run effects_LSS_file derivatives_LSS_file

import pandas as pd
from nilearn.image import load_img
from nilearn.input_data import NiftiMasker
from argparse import ArgumentParser
import numpy as np
from sklearn.covariance import ledoit_wolf
from scipy.linalg import svd


def compute_effects(args):

    # load sequence data
    allsequences = pd.read_csv(args.sequences_fn, sep = ' ')
    sequences = allsequences.loc[ allsequences.run == int(args.run), : ]
    ntypes = len(np.unique(sequences.seq_type))
    
    # open design matrix
    X_LSA = pd.read_csv(args.X_LSA_file, skiprows = 5, sep='\t', 
                        header = None).iloc[:, :-1].values

    # LSA betas
#    beta_img = load_img(args.beta_LSA_file)
#    mask = math_img('np.std(img, axis = 3) > 0', img = beta_img)
    mask_img = load_img(args.mask_file)
    masker = NiftiMasker(mask_img = mask_img)
    beta_LSA = masker.fit_transform(args.beta_LSA_file)
    res_LSA = masker.fit_transform(args.res_LSA_file)

    ntrials = sequences.shape[0]
    nconfounds = X_LSA.shape[1] - 2*ntrials

    effects_list = []
    derivatives_list = []
    for i in range(ntrials):
        print(i)
        pmat = np.zeros((X_LSA.shape[1], 2 + 2*ntypes + nconfounds))
    
        pmat[2*i, 0] = 1 # effect
        pmat[2*i + 1, 1] = 1 #derivative
        
        for j in range(ntypes):
            indices = np.where(sequences.seq_type == j + 1)[0]
            indices = indices[indices != i]
            pmat[2*indices, 2*j + 2] = 1 # effect
            pmat[2*indices + 1, 2*j + 3 ] = 1 #derivative
    
        for k in range(nconfounds):
            pmat[2*ntrials + k, 2*(ntypes + 1) + k ] = 1 # confounds
            
        A = X_LSA.transpose().dot(X_LSA)
        B = np.linalg.inv(pmat.transpose().dot(A).dot(pmat))
        C = pmat.transpose().dot(A)
        U = B.dot(C)
        
        # compute LSS effects and residuals
        beta_LSS = U.dot(beta_LSA)
        res_LSS = res_LSA + X_LSA.dot(beta_LSA - pmat.dot(beta_LSS))

        # pre-whiten
        sigma = np.std(res_LSS, axis = 0)

        if not args.multivar:
            # two first, effect and derivative
            beta_LSS_prewhitened = beta_LSS[:,:1]/sigma 
        else: 
            cov_ledoit, _ = ledoit_wolf(res_LSS)    
            Uc, Dc, Vhc = svd(cov_ledoit, full_matrices = False)
            cov_ledoit_sqrt = np.dot(Uc * np.sqrt(Dc), Vhc)
            beta_LSS_prewhitened = np.dot(beta_LSS[:,:1], 
                                          np.linalg.inv(cov_ledoit_sqrt))
    
        # we are only interested in the first parameter for each trial
        effects_list.append(beta_LSS_prewhitened[0])
        derivatives_list.append(beta_LSS_prewhitened[1])

    effects = np.vstack(effects_list)
    derivatives = np.vstack(derivatives_list)
    
    effects[np.isnan(effects)] = 0
    derivatives[np.isnan(derivatives)] = 0
    
    effects_img = masker.inverse_transform(effects)
    derivatives_img = masker.inverse_transform(derivatives)
    effects_img.to_filename(args.effects_LSS_file)
    derivatives_img.to_filename(args.derivatives_LSS_file)
    
def build_parser():
    parser = ArgumentParser()
    parser.add_argument('X_LSA_file')
    parser.add_argument('beta_LSA_file')
    parser.add_argument('res_LSA_file')
    parser.add_argument('mask_file')
    parser.add_argument('sequences_fn')
    parser.add_argument('run')
    parser.add_argument('effects_LSS_file')
    parser.add_argument('derivatives_LSS_file')
    parser.add_argument('--multivar', dest = 'multivar', 
                        action = 'store_true', 
                        required = False)
    return(parser)

def main():
    parser = build_parser()
    opts = parser.parse_args()
    compute_effects(opts)
    
if __name__== "__main__":
  main()