#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 11:58:41 2023

@author: earnestt1234
"""

import numpy as np
import pandas as pd

np.random.seed(42)

df = pd.read_csv('adni_stage_flow.csv')
df['label_from'] = df['label_from'].replace({'NS':np.nan}).astype(float)
df['label_to'] = df['label_to'].replace({'NS':np.nan}).astype(float)

label_from = df['label_from'].to_numpy()
label_to = df['label_to'].to_numpy()

def stat(label_from, label_to):
    return np.mean(label_to >= label_from)

N = 500000

observed = stat(label_from, label_to)

null = []

for n in range(N):
    swap = np.random.rand(len(df)) < 0.5
    before = np.where(swap, label_to, label_from)
    after = np.where(swap, label_from, label_to)
    null.append(stat(before, after))

null = np.array(null)
p = (null >= observed).mean()

#%% plot

import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'arial',
                     'font.size':20})

fig, ax = plt.subplots(figsize=(8, 6))

ax.hist(null, edgecolor='k', facecolor='dodgerblue', label='null',
        linewidth=3)
ax.axvline(observed, color='red', label='observed', lw=3)
ax.set_ylabel('Frequency')
ax.set_xlabel('Statistic')
ax.text(observed * 1.001, N/10, s=f'p = {round(p, 3)}', color='red', fontsize=18)
ax.axvspan(observed, null.max() * 1.010, color='red', alpha=.1)
ax.set_xlim(null.min() * .995, null.max() * 1.010)

ax.legend()

plt.tight_layout()
fig.savefig('adni_permutation_test.png', dpi=300)

