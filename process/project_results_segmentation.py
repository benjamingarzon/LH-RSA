#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 12:10:37 2021

@author: benjamin.garzon@gmail.com
"""
from nibabel.freesurfer.io import read_annot, write_morph_data, read_geometry
import pandas as pd
import numpy as np
import re, os

#annot_lh = '/data/lv0/MotorSkill/labels/GlasserParc/lh.HCP-MMP1.annot'
#annot_rh = '/data/lv0/MotorSkill/labels/GlasserParc/rh.HCP-MMP1.annot'
annot_lh = '/data/lv0/MotorSkill/parcellations/fs_LR_32/Icosahedron-162.fsaverage.L.annot'
annot_rh = '/data/lv0/MotorSkill/parcellations/fs_LR_32/Icosahedron-162.fsaverage.R.annot'

metric = 'svm_acc' #'crossnobis'
output_lh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/lh.%s.map'%metric
output_rh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/rh.%s.map'%metric
results_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/%s.csv'%metric

#results_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/crossnobis.csv'

#output_lh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/lh.alpha.map'
#output_rh = '/data/lv0/MotorSkill/fmriprep/analysis/surf/rh.alpha.map'
#results_file = '/data/lv0/MotorSkill/fmriprep/analysis/surf/alpha.csv'

def project_results_annot(annot, results_file, file_like, hemi = 'R'):
    results = pd.read_csv(results_file, sep = ',')
    print(results)
    annot_info = read_annot(annot)
    values = annot_info[0]
    output = values.astype(np.float32)*0
    labels = np.array(annot_info[2])
    results = results.loc[results.label.apply(lambda x: x[0] == hemi), :]
    for index, row in results.iterrows():
        mylabel = row.label
        if '_' in mylabel:
            mylabel = mylabel.split('_')[1]
            mylabel = re.sub('-0+', '-', mylabel)
            
        w = np.where(labels == mylabel.encode('utf-8'))[0]
        if (len(w) > 0):            
            output[values == w] = round(row.value_perm, 4)
            print(row.label, mylabel, w, round(row.value_perm, 4))
        
    write_morph_data(file_like = file_like, values = output)

project_results_annot(annot_lh, results_file, output_lh, hemi = 'L')
project_results_annot(annot_rh, results_file, output_rh, hemi = 'R')
os.system("""freeview -f """
          """/usr/local/freesurfer/7.1.1-1/subjects/fsaverage/surf/rh.inflated:"""
          """overlay=/data/lv0/MotorSkill/fmriprep/analysis/surf/rh.svm_acc.map:"""
          """annot=/data/lv0/MotorSkill/parcellations/fs_LR_32/Icosahedron-162.fsaverage.R.annot:annot_outline=1 """
          """/usr/local/freesurfer/7.1.1-1/subjects/fsaverage/surf/lh.inflated:"""
          """overlay=/data/lv0/MotorSkill/fmriprep/analysis/surf/lh.svm_acc.map:"""
          """annot=/data/lv0/MotorSkill/parcellations/fs_LR_32/Icosahedron-162.fsaverage.L.annot:annot_outline=1""")

#freeview -f /usr/local/freesurfer/7.1.1-1/subjects/fsaverage/surf/lh.inflated:overlay=/data/lv0/MotorSkill/labels/GlasserParc/lh.test:annot=lh.HCP-MMP1.annot 
