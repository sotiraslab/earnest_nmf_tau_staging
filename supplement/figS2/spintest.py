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
from neuromaps.nulls import alexander_bloch, vasa
from neuromaps.stats import compare_images
import nibabel as nib
import pandas as pd
from scipy.io import loadmat
from statsmodels.stats.multitest import multipletests

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

rotated = alexander_bloch(df['ADNI'],
                          atlas='fsaverage',
                          density='10k',
                          n_perm=N_PERM,
                          seed=42,
                          parcellation=parcellation)

def metric(a, b):
    return 1.0

from scipy.spatial.distance import cosine

def myfunc(a, b):
    return np.float64(cosine(a, b))

corr, pval = compare_images(df['ADNI'],
                            df['OASIS'],
                            metric = myfunc,
                            nulls=rotated)
