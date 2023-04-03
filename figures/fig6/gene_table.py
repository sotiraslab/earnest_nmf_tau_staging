#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 09:38:57 2023

@author: earnestt1234
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('wightman_hits.csv')
table = df.pivot(index='Gene', columns='Component', values='R')
table = table.reindex(['MedialTemporal', 'Occipital', 'Sensorimotor', 'Orbitofrontal'], axis=1)
rows, cols = table.shape

plt.rcParams.update({'font.family': 'arial'})

plt.figure(figsize=(4.5, 10), dpi=300)
plt.imshow(table, cmap='bwr', vmin=-1, vmax=1, aspect='auto')
plt.grid(which='minor')

ax = plt.gca()

ax.xaxis.tick_top()
ax.set_xticks(np.arange(table.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(table.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="gray", linestyle='-')
ax.tick_params(which="both", bottom=False, left=False, top=False)

for spine in ax.spines.values():
        spine.set_edgecolor('gray')

plt.yticks(range(rows), table.index, weight='bold')
plt.xticks(range(cols), ['PTC1-MedialTemporal',
                         'PTC5-Occipital',
                         'PTC7-Sensorimotor',
                         'PTC8-Orbitofrontal'],
           ha='left', rotation=30, weight='bold')

for i in range(rows):
    for j in range(cols):
        val = table.iloc[i, j]
        if pd.isna(val):
            continue
        plt.text(j, i, '{0:.3f}'.format(val), ha='center', va='center', fontsize=12, color='white')

# save
plt.savefig('gene_table.png', bbox_inches='tight')

# ----------
# Colorbar
# ----------

# gradient = np.linspace(0, 1, 256)
# gradient = np.vstack((gradient, gradient))

# fig, ax = plt.subplots(figsize=(7, 1), dpi=300)
# ax.imshow(gradient, aspect='auto', cmap='bwr')
# ax.axis('off')

# fig.savefig('bwr_colorbar.png', transparent=True)
