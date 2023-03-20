#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 09:34:46 2023

@author: earnestt1234
"""

import pandas as pd

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({'font.family':'arial',
                     'font.size': 20})

# read
df = pd.read_csv('adni_stage_flow.csv')
labels = sorted(df['label_to'].unique(), key = lambda x: (not x.isnumeric(), x))
# labels = [l for l in labels if l.isnumeric()]

# create heatmap data
transitions = df.groupby(['label_from', 'label_to'])['subject'].count()
new_index = pd.MultiIndex.from_product([labels, labels])
transitions = transitions.reindex(new_index).reset_index()
transitions.columns = ['label_from', 'label_to', 'n']
transitions['n'] = transitions['n'].fillna(0)
transitions['perc'] = transitions.groupby('label_from')['n'].transform(lambda x: (x / x.sum()) * 100).fillna(0)

tbl = transitions.pivot(index='label_from', columns='label_to', values='perc')

gs = plt.GridSpec(nrows=1, ncols=2, width_ratios=[20, 1], wspace=0)

fig = plt.figure(figsize=(9, 8))
ax = fig.add_subplot(gs[0])
c1 = fig.add_subplot(gs[1])

ax.set_aspect('equal')
sns.heatmap(tbl, ax=ax, cmap='inferno', vmin=0, vmax=100, cbar_ax=c1)

ax.tick_params(axis='both', which='both',length=0)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.set_ylabel('Stage (scan 1)')
ax.set_xlabel('Stage (scan 2)')
ax.set_xlim([-.1, 6.1])
ax.set_ylim([6.1, -.1])

# draw boxes on diagonal
for i in range(len(labels)):
    for j in range(len(labels)):
        if i == j:
            rec = patches.Rectangle((i, j), width=1, height=1, color='darkgray', fill=False, lw=5, zorder=5)
            ax.add_patch(rec)

# draw NS box
r1 = patches.Rectangle((0, 5), width=6, height=1, color='white', fill=False, lw=5, zorder=1)
r2 = patches.Rectangle((5, 0), width=1, height=6, color='white', fill=False, lw=5, zorder=1)
ax.add_patch(r1)
ax.add_patch(r2)

plt.savefig('adni_longitudinal_stage_heatmap.png', dpi=300)
