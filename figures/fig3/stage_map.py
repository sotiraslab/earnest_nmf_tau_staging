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

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----------
# required files
# ----------

PATH_WTA = '../../derivatives/adni/ptc8_winner_take_all.csv'
PATH_SCRIPTS = '../../scripts'

# ----------
# setup
# ----------

stages = {'Cmp.MedialTemporal': 1,
          'Cmp.LeftParietalTemporal': 2,
          'Cmp.RightParietalTemporal': 2,
          'Cmp.Precuneus': 2,
          'Cmp.Occipital': 3,
          'Cmp.LateralFrontal': 3,
          'Cmp.Sensorimotor': 4,
          'Cmp.Orbitofrontal': 4}

dkt_table = pd.read_csv(PATH_WTA)
dkt_table['region'] = dkt_table['label']
dkt_table['value'] = dkt_table['name'].map(stages)

# create colormap
colors = np.zeros((256, 4))
colors[:65, :] = [0, 158, 115, 255]
colors[65:129, :] = [240, 226, 66, 255]
colors[129:193, :] = [230, 159, 0, 255]
colors[193:, :] = [213, 94, 0, 255]
colors /= 255
cmap = mpl.colors.ListedColormap(colors, name='stages')

try:
    _ = mpl.colormaps.get_cmap('stages')
except ValueError:
    mpl.colormaps.register(cmap)

# ----------
# read plotting func
# ----------

sys.path.append(PATH_SCRIPTS)

from plot_brainsapce import plot_dkt_table_brainspace

#%%
# ----------
# plot all stages
# ----------


plot_dkt_table_brainspace(dkt_table,
                          layer='pial',
                          size=(800, 500),
                          cmap='stages',
                          nan_color=(0.5, 0.5, 0.5, 1),
                          layout_style='grid',
                          filename='stage_map.png',
                          screenshot=True,
                          zoom=1.8)

dkt_table.to_csv('../../derivatives/adni/ptc8_stage_assignment.csv')

#%%
# ----------
# plot individual stages
# ----------

stage_colors = [[0, 158, 115, 255],
                [240, 226, 66, 255],
                [230, 159, 0, 255],
                [213, 94, 0, 255]]

for i in [1, 2, 3, 4]:
    tmp = dkt_table.copy()
    tmp.loc[tmp['value'] != i, 'value'] = None

    tmp_map = np.zeros((256, 4))
    tmp_map[:, :] = stage_colors[i-1]
    tmp_map /= 255
    tmp_map = mpl.colors.ListedColormap(tmp_map, name=f'stage_{i}')

    try:
        _ = mpl.colormaps.get_cmap(f'stage_{i}')
    except ValueError:
        mpl.colormaps.register(tmp_map)

    plot_dkt_table_brainspace(tmp,
                              layer='pial',
                              size=(800, 500),
                              cmap=f'stage_{i}',
                              nan_color=(1.0, 1.0, 1.0, 1),
                              layout_style='grid',
                              filename=f'stage_map_{i}.png',
                              screenshot=True,
                              zoom=1.8)



#%%
# ----------
# colormap
# ----------

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

fig, ax = plt.subplots(figsize=(7, 1), dpi=300)
ax.imshow(gradient, aspect='auto', cmap='stages')
ax.axis('off')

for i in range(4):
    ax.text(x=64*i + 32, y=0.5, s=i+1, ha='center', va='center', font='arial', size=20)

fig.savefig('colorbar.png', transparent=True)

