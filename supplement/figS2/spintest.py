#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 15:11:22 2023

@author: tom.earnest
"""

# ----------
# imports
# ----------

import abagen
import numpy as np
from neuromaps.nulls import vasa
from neuromaps.stats import compare_images
import nibabel as nib
import pandas as pd
from scipy.special import comb
from scipy.io import loadmat

# ----------
# required files
# ----------

# must be run after figure 6!
FIG6_TEMPLATE = '../../figures/fig6/ptcs_prepped.csv'

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
adni_w = loadmat(PATH_NMF_MAT_ADNI)['Wnorm']
oasis_w = loadmat(PATH_NMF_MAT_OASIS)['Wnorm']
df = pd.read_csv(FIG6_TEMPLATE)

# ---------
# update table
# ---------

adni = pd.DataFrame({'hemisphere': regions.str.extract('(lh|rh)', expand=False),
                     'label': regions.str.replace('(lh|rh)_', '', regex=True),
                     'ADNI': adni_w.argmax(axis=1).astype(float)})
adni['hemisphere'] = adni['hemisphere'].map({'lh':'L', 'rh': 'R'})

oasis = pd.DataFrame({'hemisphere': regions.str.extract('(lh|rh)', expand=False),
                     'label': regions.str.replace('(lh|rh)_', '', regex=True),
                     'OASIS': oasis_w.argmax(axis=1).astype(float)})
oasis['hemisphere'] = oasis['hemisphere'].map({'lh':'L', 'rh': 'R'})
matches = loadmat(PATH_PTC_MATCH)['idx_hug1']
mapper = dict(zip(matches.astype(float)[:, 0] - 1, np.arange(K).astype(float)))
oasis['OASIS'] = oasis['OASIS'].map(mapper)

df = df.merge(adni, how='left', on = ['hemisphere', 'label'])
df = df.merge(oasis, how='left', on = ['hemisphere', 'label'])

# -------
# prep neuromaps
# --------

# dkt atlas
dkt = abagen.fetch_desikan_killiany(surface=True)
dkt_l_path, dkt_r_path = dkt['image']
dkt_info_path = dkt['info']

dkt_l = nib.load(dkt_l_path)
dkt_r =  nib.load(dkt_r_path)
parcellation = (dkt_l, dkt_r)

dkt_labels = pd.read_csv(dkt_info_path)
dkt_labels = dkt_labels[dkt_labels['structure'] == 'cortex']

# ---------
# run
# ---------

N_PERM = 10000

rotated = vasa(df['ADNI'],
               atlas='fsaverage',
               density='10k',
               n_perm=N_PERM,
               seed=42,
               parcellation=parcellation)

def adjusted_rand_index(x, y):

    n = len(x)

    clusters_x = np.unique(x)
    clusters_y = np.unique(y)
    cmat = np.zeros((len(clusters_x), len(clusters_y)))

    for i, xlabel in enumerate(clusters_x):
        for j, ylabel in enumerate(clusters_y):
            cmat[i, j] = np.logical_and(x == xlabel, y == ylabel).sum()

    a = cmat.sum(axis=1)
    b = cmat.sum(axis=0)

    cmat_comb = comb(cmat, 2).sum()
    a_comb = comb(a, 2).sum()
    b_comb = comb(b, 2).sum()
    n_comb = comb(n, 2)

    top = cmat_comb - ((a_comb * b_comb)/n_comb)
    bot = (0.5 * (a_comb + b_comb)) - ((a_comb * b_comb)/n_comb)

    ari = top/bot
    return np.float64(ari)

ari, pval, nulls  = compare_images(df['ADNI'],
                                   df['OASIS'],
                                   metric = adjusted_rand_index,
                                   ignore_zero=False,
                                   nulls=rotated,
                                   return_nulls=True)

#%% save

string = f'observed ARI = {ari}, p-value = {pval}'
with open('spintest_results.txt', 'w') as f:
    f.write(string)

#%% plot

import matplotlib.pyplot as plt

plt.rcParams.update({'font.family':'arial',
                     'font.size':20})

fig, ax = plt.subplots(figsize=(8, 6))

ax.hist(nulls, edgecolor='k', facecolor='dodgerblue', label='null',
        linewidth=3)
ax.axvline(ari, color='red', label='observed', lw=3)
ax.set_ylabel('Frequency')
ax.set_xlabel('Adjusted Rand Index')
# ax.text(ari * 1.001, N_PERM/10, s=f'p = {round(pval, 3)}', color='red', fontsize=18)
ax.set_xlim(nulls.min() * .995, nulls.max() * 1.1)

ax.legend()

plt.tight_layout()
fig.savefig('winner_take_all_spintest.png', dpi=300)
