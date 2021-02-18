#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 12:10:37 2021

@author: benjamin.garzon@gmail.com
"""
from nibabel.freesurfer.io import read_annot, write_morph_data, read_geometry
import pandas as pd
import numpy as np

annot_lh = '/data/lv0/MotorSkill/labels/GlasserParc/lh.HCP-MMP1.annot'
annot_rh = '/data/lv0/MotorSkill/labels/GlasserParc/rh.HCP-MMP1.annot'

output_lh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/lh.svm.map'
output_rh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/rh.svm.map'
results_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/svm_acc.csv'

output_lh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/lh.alpha.map'
output_rh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/rh.alpha.map'
results_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/alpha.csv'

def project_results_annot(annot, results_file, file_like):
    results = pd.read_csv(results_file, sep = ',')
    print(results)
    annot_info = read_annot(annot)
    values = annot_info[0]
    output = values.astype(np.float32)*0
    labels = np.array(annot_info[2])
    for index, row in results.iterrows():
        w = np.where(labels == row.label.encode('utf-8'))[0]
        if (len(w) > 0):            
            output[values == w] = row.value
            print(row.label, w, row.value)
        
    write_morph_data(file_like = file_like, values = output)

project_results_annot(annot_lh, results_file, output_lh)
project_results_annot(annot_rh, results_file, output_rh)
#freeview -f /usr/local/freesurfer/7.1.1-1/subjects/fsaverage/surf/lh.inflated:overlay=/data/lv0/MotorSkill/labels/GlasserParc/lh.test:annot=lh.HCP-MMP1.annot 
