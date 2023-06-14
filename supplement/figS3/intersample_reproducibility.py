#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 12:35:08 2022

@author: earnestt1234
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

#%% Load comparison results

path = 'adni_v_oasis_compare/reproducibilityResults.mat'

d = loadmat(path)
ari = d['ARI']
meanInner = d['meanInner']
meanInner = np.concatenate(*meanInner)
medianInner = d['medianInner']
medianInner = np.concatenate(*medianInner)
ranks = d['sortedBasisNum'].flatten()
xticks = [i for i in ranks if i % 2 == 0]


#%% Plot params

plt.rcParams.update({'font.family':'arial',
                     'font.size': 24})

#%% ari mean
x = ranks

color = 'dodgerblue'
ylab = 'Adjusted Rand Index'
xlab = 'Number of PTCs'

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(x, ari, zorder=3, lw=3, label='ADNI-ADS/OASIS3-ADS', color='darkviolet')


ax.set_xticks(xticks)
ax.grid(zorder=1, alpha=.3)
ax.set_ylabel(ylab)
ax.set_xlabel(xlab)
ax.legend()

plt.tight_layout()
plt.savefig('matched_ari.png', dpi=300)

#%% mean inner

x = ranks
ylab = 'Mean Inner Product'
xlab = 'Number of PTCs'

fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(x, meanInner, zorder=3, lw=3, label='ADNI-ADS/OASIS3-ADS', color='darkviolet')
ax.set_xticks(xticks)
ax.grid(zorder=1, alpha=.3)
ax.set_ylabel(ylab)
ax.set_xlabel(xlab)
ax.legend()

plt.tight_layout()
plt.savefig('matched_mean.png', dpi=300)
