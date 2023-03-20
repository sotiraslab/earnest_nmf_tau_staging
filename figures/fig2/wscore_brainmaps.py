#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:05:59 2023

@author: tom.earnest
"""

# ----------
# imports
# ----------

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----------
# required files
# ----------

PATH_WTA = '../../derivatives/adni/ptc8_winner_take_all.csv'
PATH_ADNI_WORDER = '../../derivatives/adni/wscore_stage_order.csv'
PATH_ADNI_WSCORES = '../../derivatives/adni/data_with_wscores.csv'
PATH_OASIS_WORDER = '../../derivatives/oasis3/wscore_stage_order.csv'
PATH_OASIS_WSCORES = '../../derivatives/oasis3/data_with_wscores.csv'
PATH_SCRIPTS = '../../scripts'

# ----------
# read plotting func
# ----------

sys.path.append(PATH_SCRIPTS)

from plot_brainsapce import plot_dkt_table_brainspace


# ----------
# ADNI
# ----------

wta = pd.read_csv(PATH_WTA)
worder = pd.read_csv(PATH_ADNI_WORDER)
wscores = pd.read_csv(PATH_ADNI_WSCORES)
wscores = wscores[wscores['Group'] == 'TrainingBaseline']
wdf = wscores[wscores.columns[wscores.columns.str.contains('WScore')]]
n_any_elevated = (wdf >= 2.5).any(axis=1).sum()
worder['NPos'] = (worder['NPos'] / n_any_elevated) * 100
worder['Region'] = worder['Region'].str.replace('.WScore', '', regex=False)

dkt_table = wta.merge(worder, how='left', left_on='name', right_on='Region')
dkt_table['region'] = dkt_table['label']
dkt_table['value'] = dkt_table['NPos']

fig = plot_dkt_table_brainspace(dkt_table,
                                layer='pial',
                                size=(800, 500),
                                cmap='Reds',
                                color_range=(20, 90),
                                nan_color=(0.5, 0.5, 0.5, 1),
                                layout_style='grid',
                                filename='adni_w_map.png',
                                screenshot=True,
                                zoom=1.8)

# ----------
# OASIS
# ----------

wta = pd.read_csv(PATH_WTA)
worder = pd.read_csv(PATH_OASIS_WORDER)
wscores = pd.read_csv(PATH_OASIS_WSCORES)
wscores = wscores[wscores['Group'] == 'TrainingSet']
wdf = wscores[wscores.columns[wscores.columns.str.contains('WScore')]]
n_any_elevated = (wdf >= 2.5).any(axis=1).sum()
worder['NPos'] = (worder['NPos'] / n_any_elevated) * 100
worder['Region'] = worder['Region'].str.replace('.WScore', '', regex=False)

dkt_table = wta.merge(worder, how='left', left_on='name', right_on='Region')
dkt_table['region'] = dkt_table['label']
dkt_table['value'] = dkt_table['NPos']

fig = plot_dkt_table_brainspace(dkt_table,
                                layer='pial',
                                size=(800, 500),
                                cmap='Reds',
                                color_range=(20, 90),
                                nan_color=(0.5, 0.5, 0.5, 1),
                                layout_style='grid',
                                filename='oasis_w_map.png',
                                screenshot=True,
                                zoom=1.8)


# ----------
# colormap
# ----------

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

fig, ax = plt.subplots(figsize=(7, 1), dpi=300)
ax.imshow(gradient, aspect='auto', cmap='Reds')
ax.axis('off')

fig.savefig('colorbar.png', transparent=True)
