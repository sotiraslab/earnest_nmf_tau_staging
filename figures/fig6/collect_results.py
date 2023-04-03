#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 10:16:00 2022

@author: tom.earnest
"""

import os
import pandas as pd

wightman_genes = []

csvs = sorted(os.listdir(os.getcwd()))
csvs = [csv for csv in csvs if 'genecorrelation' in csv]
names = [s.replace('genecorrelation_', '').replace('.csv', '') for s in csvs]

for i, csv in enumerate(csvs):

    df = pd.read_csv(csv)
    hits = df['fdr_bh'].le(0.05).sum()

    hits_df = df[df['fdr_bh'].le(0.05)]
    for j, row in hits_df.iterrows():
        x = {'Component':names[i],
             'Gene':row['gene'],
             'R': round(row['corr'], 3),
             'p': round(row['fdr_bh'], 3)}
        wightman_genes.append(x)

wightman_genes = pd.DataFrame(wightman_genes).sort_values(['Component', 'Gene'])
wightman_genes['Component'] = wightman_genes['Component'].replace(
    {'Cmp.MedialTemporal': 'PTC1-MedialTemporal',
     'Cmp.Occipital': 'PTC5-Occipital',
     'Cmp.OrbitoFrontal': 'PTC8-Orbitofrontal',
     'Cmp.SensoryMotor': 'PTC7-Sensorimotor'})

wightman_genes.to_csv('wightman_hits.csv', index=False)
