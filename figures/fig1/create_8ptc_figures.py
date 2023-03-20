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

# ----------
# required files
# ----------

PATH_REGIONS = '../../derivatives/adni/nmf_regions_ggseg.csv'
PATH_PTC_NAMES = '../../derivatives/adni/names_8ptc.csv'
PATH_NMF_MAT = '../../nmf/adni/results/mat/NumBases8.mat'
PATH_SCRIPTS = '../../scripts'

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
# plot
# ----------


for i, label in enumerate(ptc_names):
    print(f'({i}) plotting PTC-{label}...')
    dkt_table = nmf_component_to_dkt_table(PATH_NMF_MAT, component_index=i, region_names=regions)
    plot_dkt_table_brainspace(dkt_table,
                              layer='pial',
                              size=(1600, 300),
                              cmap='plasma',
                              nan_color=(0.5, 0.5, 0.5, 1),
                              filename=f'{label}.png',
                              screenshot=True,
                              zoom=1.7)
