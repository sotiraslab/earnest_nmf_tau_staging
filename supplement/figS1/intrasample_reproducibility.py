#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 12:35:08 2022

@author: earnestt1234
"""

import matplotlib.pyplot as plt
from scipy.io import loadmat

#%% Load ADNI

path = '../../nmf/adni/results/repeat_reproducibility_mat/reproducibilityResults.mat'

d = loadmat(path)
ari_ADNI = d['repeat_ARI']
meanInner_ADNI = d['repeat_meanInner']
medianInner_ADNI = d['repeat_medianInner']
ranks = d['sortedBasisNum'].flatten()
xticks = [i for i in ranks if i % 2 == 0]

#%% Load ADRC

path = '../../nmf/oasis3/results/repeat_reproducibility_mat/reproducibilityResults.mat'

d = loadmat(path)
ari_ADRC = d['repeat_ARI']
meanInner_ADRC = d['repeat_meanInner']
medianInner_ADRC = d['repeat_medianInner']
ranks = d['sortedBasisNum'].flatten()
xticks = [i for i in ranks if i % 2 == 0]


#%% Plot params

plt.rcParams.update({'font.family':'arial',
                     'font.size': 24})

#%% ari mean
x = ranks

err = ari_ADNI.std(axis=0)
color = 'dodgerblue'
ylab = 'Adjusted Rand Index'
xlab = 'Number of PTCs'

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(x, ari_ADNI.mean(axis=0), zorder=3, lw=3, label='ADNI-ADS', color='red')
ax.plot(x, ari_ADRC.mean(axis=0), zorder=3, lw=3, label='OASIS3-ADS', color='blue', linestyle='dashed')

ax.set_xticks(xticks)
ax.grid(zorder=1, alpha=.3)
ax.set_ylabel(ylab)
ax.set_xlabel(xlab)
ax.legend()

plt.tight_layout()
plt.savefig('ari_mean.png', dpi=300)

#%% mean inner

x = ranks
ylab = 'Mean Inner Product'
xlab = 'Number of PTCs'

fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(x, meanInner_ADNI.mean(axis=0), zorder=3, lw=3, label='ADNI-ADS', color='red')
ax.plot(x, meanInner_ADRC.mean(axis=0), zorder=3, lw=3, label='OASIS3-ADS', color='blue', linestyle='dashed')
ax.set_xticks(xticks)
ax.grid(zorder=1, alpha=.3)
ax.set_ylabel(ylab)
ax.set_xlabel(xlab)
ax.legend()

plt.tight_layout()
plt.savefig('mean_inner_product.png', dpi=300)
