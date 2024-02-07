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

import matplotlib as mpl
import numpy as np
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

from plot_brainsapce import nmf_winner_take_all_dkt_table, plot_dkt_table_brainspace

# ----------
# colormap
# ----------

def get_eighth(n):
    return np.arange(n*32, (n*32) + 32)

# create colormap
colors = np.zeros((256, 4))
colors[get_eighth(0), :] = [255,  58,  58, 255]
colors[get_eighth(1), :] = [252, 149,  59, 255]
colors[get_eighth(2), :] = [255, 232,  30, 255]
colors[get_eighth(3), :] = [ 73, 212,  64, 255]
colors[get_eighth(4), :] = [ 52, 236, 230, 255]
colors[get_eighth(5), :] = [ 14, 113, 189, 255]
colors[get_eighth(6), :] = [145,  35, 209, 255]
colors[get_eighth(7), :] = [142,  96,  61, 255]
colors /= 255
cmap = mpl.colors.ListedColormap(colors, name='stages')

if 'stages' in list(mpl.colormaps.keys()):
    mpl.colormaps.unregister('stages')
    
mpl.colormaps.register(cmap)

# ----------
# plot ADNI
# ----------

dkt_table = nmf_winner_take_all_dkt_table(PATH_NMF_MAT_ADNI,
                                          region_names=regions)
plot_dkt_table_brainspace(dkt_table,
                          layer='pial',
                          size=(1600, 300),
                          cmap='stages',
                          nan_color=(0.5, 0.5, 0.5, 1),
                          filename='WTA_ADNI.png',
                          screenshot=True,
                          zoom=1.7)

# ----------
# plot OASIS
# ----------


dkt_table = nmf_winner_take_all_dkt_table(PATH_NMF_MAT_OASIS,
                                          region_names=regions)

# make ordering of components match adni, based on hungarian
matches = loadmat(PATH_PTC_MATCH)['idx_hug1']
mapper = dict(zip(matches.astype(float)[:, 0] - 1, np.arange(K)))
dkt_table['value'] = dkt_table['value'].map(mapper)
oasis = dkt_table.copy()
plot_dkt_table_brainspace(dkt_table,
                          layer='pial',
                          size=(1600, 300),
                          cmap='stages',
                          nan_color=(0.5, 0.5, 0.5, 1),
                          filename='WTA_OASIS.png',
                          screenshot=True,
                          zoom=1.7)
