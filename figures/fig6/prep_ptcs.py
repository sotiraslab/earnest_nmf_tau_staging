#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:36:56 2022

@author: tom.earnest
"""

#%% imports

import os

import pandas as pd
from scipy.io import loadmat

#%% Required files

PATH_8PTC = '../../nmf/adni/results/mat/NumBases8.mat'
PATH_8PTC_NAMES = '../../derivatives/adni/names_8ptc.csv'
PATH_DKT_LABELS = 'dkt_labels.csv'
PATH_NMF_REGIONS = '../../derivatives/adni/nmf_regions.csv'

#%% read components

path = os.path.abspath(PATH_8PTC)
mat = loadmat(path)
wnorm = mat['Wnorm']

wnames = pd.read_csv(PATH_8PTC_NAMES)['Component'].to_list()

#%% read DKT labels

dkt_labels = pd.read_csv('dkt_labels.csv')
dkt_labels = dkt_labels[dkt_labels['structure'] == 'cortex']

#%% Order regions in the abagen parcellation order

regions = pd.read_csv(PATH_NMF_REGIONS)['Feature']
regions = regions.str.lower().str.replace('ctx_','').str.replace('_suvr', '').str.split('_', expand=True)
regions.columns = ['hemisphere','label']

# recode some things to map with abagen
regions['hemisphere'] = regions['hemisphere'].map({'lh':'L', 'rh':'R'})

wdf = pd.concat([regions, pd.DataFrame(wnorm, columns=wnames)], axis=1)
wdf_ordered = dkt_labels.merge(wdf, how='outer', on=['hemisphere', 'label'])
wdf_ordered.to_csv('ptcs_prepped.csv', index=False)

