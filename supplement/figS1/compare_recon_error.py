#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 16:13:42 2022

@author: earnestt1234
"""

import os

import matplotlib.pyplot as plt
import natsort
import numpy as np
from scipy.io import loadmat

#%% Load ADNI

X_path = '../../derivatives/adni/nmf_matrix.csv'
matfolder = '../../nmf/adni/results/mat'

X = np.loadtxt(X_path, delimiter=',')

errors = []
for mat in natsort.natsorted(os.listdir(matfolder)):
    if not mat.lower().endswith('mat'): continue;
    print(f'{mat}...')
    fullfile = os.path.join(matfolder, mat)
    d = loadmat(fullfile)
    W = d['W']
    H = d['H']
    fro = np.linalg.norm(X - np.matmul(W, H), ord='fro')
    errors.append(fro)

adni_errors = np.array(errors)
ranks = np.array(range(2, 21))

#%% Load ADRC

X_path = '../../derivatives/oasis3/nmf_matrix.csv'
matfolder = '../../nmf/oasis3/results/mat'

X = np.loadtxt(X_path, delimiter=',')

errors = []
for mat in natsort.natsorted(os.listdir(matfolder)):
    if not mat.lower().endswith('mat'): continue;
    print(f'{mat}...')
    fullfile = os.path.join(matfolder, mat)
    d = loadmat(fullfile)
    W = d['W']
    H = d['H']
    fro = np.linalg.norm(X - np.matmul(W, H), ord='fro')
    errors.append(fro)

adrc_errors = np.array(errors)
ranks = np.array(range(2, 21))

#%% Normalize errors

adni_errors_norm = (adni_errors - adni_errors.min()) / (adni_errors.max() - adni_errors.min())
adrc_errors_norm = (adrc_errors - adrc_errors.min()) / (adrc_errors.max() - adrc_errors.min())

#%% 1. Recon error

plt.rcParams.update({'font.family':'arial',
                     'font.size': 24})
xticks = range(2, 21, 2)

x = ranks

ylab = 'Scaled reconstruction error'
xlab = 'Number of PTCs'

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(ranks, adni_errors_norm, zorder=3, color='red', lw=3, label='ADNI-ADS')
ax.plot(ranks, adrc_errors_norm, zorder=3, color='blue', lw=3, label='OASIS3-ADS', linestyle='dashed')

ax.set_xticks(xticks)
ax.grid(zorder=1, alpha=.3)
ax.set_ylabel(ylab)
ax.set_xlabel(xlab)
ax.legend()

plt.tight_layout()
plt.savefig('compare_recon_error.png', dpi=300)

#%% 2. Gradient recon error

# this is using the surrounding points, instead of the
# difference with previous

x = ranks
y = np.gradient(errors)
ylab = 'Gradient of scaled recon. error'
xlab = 'Number of PTCs'

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(ranks, np.gradient(adni_errors_norm), zorder=3, color='red', lw=3, label='ADNI-ADS')
ax.plot(ranks, np.gradient(adrc_errors_norm), zorder=3, color='blue', lw=3, label='OASIS3-ADS', linestyle='dashed')

ax.set_xticks(xticks)
ax.grid(zorder=1, alpha=.3)
ax.set_ylabel(ylab)
ax.set_xlabel(xlab)
ax.legend()

plt.tight_layout()
fig.savefig('compare_gradient.png', dpi=300)
