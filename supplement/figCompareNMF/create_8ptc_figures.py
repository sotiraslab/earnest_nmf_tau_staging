#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:09:12 2023

@author: earnestt1234
"""

# ----------
# imports
# ----------

import sys

import pandas as pd
from scipy.io import loadmat

# ----------
# required files
# ----------

K = 8

PATH_REGIONS = '../../derivatives/adni/nmf_regions_ggseg.csv'
PATH_PTC_NAMES = f'../../derivatives/adni/names_{K}ptc.csv'
PATH_NMF_MAT_ADNI = f'../../nmf/adni/results/mat/NumBases{K}.mat'
PATH_NMF_MAT_OASIS = f'../../nmf/oasis3/results/mat/NumBases{K}.mat'
PATH_SCRIPTS = '../../scripts'
PATH_PTC_MATCH = f'adni_v_oasis_matching/matching/Match{K}.mat'

# ----------
# read inputs
# ----------

regions = pd.read_csv(PATH_REGIONS)['label']
ptc_names = pd.read_csv(PATH_PTC_NAMES)['Component']

# ----------
# read plotting func
# ----------

sys.path.append(PATH_SCRIPTS)

from plot_brainsapce import nmf_component_to_dkt_table, plot_dkt_table_brainspace

# ----------
# plot ADNI
# ----------


for i, label in enumerate(ptc_names):
    print(f'({i}) plotting PTC-{label}...')
    dkt_table = nmf_component_to_dkt_table(PATH_NMF_MAT_ADNI, component_index=i, region_names=regions)
    plot_dkt_table_brainspace(dkt_table,
                              layer='pial',
                              size=(1600, 300),
                              cmap='plasma',
                              nan_color=(0.5, 0.5, 0.5, 1),
                              filename=f'{label}-ADNI.png',
                              screenshot=True,
                              zoom=1.7)

# ----------
# plot OASIS
# ----------

matches = loadmat(PATH_PTC_MATCH)['idx_hug1']

for i, label in enumerate(ptc_names):
    print(f'({i}) plotting OASIS match for PTC-{label}...')
    component_index = matches[i][0] - 1
    dkt_table = nmf_component_to_dkt_table(PATH_NMF_MAT_OASIS, component_index=component_index, region_names=regions)
    plot_dkt_table_brainspace(dkt_table,
                              layer='pial',
                              size=(1600, 300),
                              cmap='plasma',
                              nan_color=(0.5, 0.5, 0.5, 1),
                              filename=f'{label}-OASIS.png',
                              screenshot=True,
                              zoom=1.7)
