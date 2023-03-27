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

PATH_REGIONS = '../../derivatives/adni/nmf_regions_ggseg.csv'
PATH_NMF_MAT_ADNI = '../../nmf/adni/results/mat/NumBases2.mat'
PATH_NMF_MAT_OASIS = '../../nmf/oasis3/results/mat/NumBases2.mat'
PATH_SCRIPTS = '../../scripts'

# matched components are provided in the figS1 directory
PATH_2PTC_MATCH = '../../supplement/figS1/adni_v_oasis_compare/matching/Match2.mat'

# ----------
# read inputs
# ----------

regions = pd.read_csv(PATH_REGIONS)['label']

# ----------
# read plotting func
# ----------

sys.path.append(PATH_SCRIPTS)

from plot_brainsapce import nmf_component_to_dkt_table, plot_dkt_table_brainspace

# ----------
# plot ADNI
# ----------


for i in range(2):
    print(f'({i}) plotting PTC-{i+1}...')
    dkt_table = nmf_component_to_dkt_table(PATH_NMF_MAT_ADNI, component_index=i, region_names=regions)
    plot_dkt_table_brainspace(dkt_table,
                              size=(800, 500),
                              layer='pial',
                              cmap='plasma',
                              nan_color=(0.5, 0.5, 0.5, 1),
                              layout_style='grid',
                              filename=f'ptc{i+1}-ADNI.png',
                              screenshot=True,
                              zoom=1.8)

# ----------
# plot OASIS
# ----------

matches = loadmat(PATH_2PTC_MATCH)['idx_hug1']

for i in range(2):
    print(f'({i}) plotting OASIS match for PTC-{i+1}...')
    component_index = matches[i][0] - 1
    dkt_table = nmf_component_to_dkt_table(PATH_NMF_MAT_OASIS, component_index=component_index, region_names=regions)
    plot_dkt_table_brainspace(dkt_table,
                              layer='pial',
                              size=(800, 500),
                              cmap='plasma',
                              nan_color=(0.5, 0.5, 0.5, 1),
                              filename=f'ptc{i+1}-OASIS.png',
                              layout_style='grid',
                              screenshot=True,
                              zoom=1.8)
