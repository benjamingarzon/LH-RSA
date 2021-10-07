#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 14:56:24 2021

@author: xgarzb@GU.GU.SE
"""
import os
import numpy as np
from roi_analysis_funcs import get_across_session_scores, process_scores
analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis'
suffix = 'mask-cross-runprew'
permutate = False ## fix permutation method
n_sample = 10
output_path = os.path.join(analysis_dir, 
                                'surf', 
                                'across_session_scores_complete_%s.csv'%(suffix))

scores_df = get_across_session_scores(os.path.join(analysis_dir, 
                                'roi_data'),
                                      suffix,
                                      output_path = output_path,
                                      num_cores = 35,
#                                      subject_mask_remove = 'sub-lue1',
#                                      selected_subjects = ['sub-lue2103', 'sub-lue2203',
#                                                           'sub-lue3103', 'sub-lue3203',
#                                                           'sub-lue5103', 'sub-lue5203',
#                                                           'sub-lue4103', 'sub-lue4203'],
#                                      selected_sessions = [1, 2, 3, 4],
#                                      selected_labels = ['R_SPL', 'R_C1'],
                                      permutate = permutate,
                                      n_sample = n_sample)

