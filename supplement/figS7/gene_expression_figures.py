#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:09:12 2023

@author: earnestt1234
"""

# ----------
# imports
# ----------

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----------
# required files
# ----------

# THIS MUST BE RUN AFTER TABLE 2 HAS BEEN RUN

PATH_ABAGEN_EXPRESSION = '../../tables/table2/abagen_expression_dkt.csv'
PATH_DKT_LABELS = '../../tables/table2/dkt_labels.csv'
PATH_HITS = '../../tables/table2/wightman_hits.csv'
PATH_SCRIPTS = '../../scripts'

# ----------
# read inputs
# ----------

expression = pd.read_csv(PATH_ABAGEN_EXPRESSION)
dkt_labels = pd.read_csv(PATH_DKT_LABELS)
wightman_hits = pd.read_csv(PATH_HITS)

expression['id'] = expression['label']

dkt_labels = dkt_labels.loc[dkt_labels['structure'] == 'cortex', :]
dkt_labels['region'] = dkt_labels['hemisphere'].str.lower() + 'h_' + dkt_labels['label']

# ----------
# read plotting func
# ----------

sys.path.append(PATH_SCRIPTS)

from plot_brainsapce import plot_dkt_table_brainspace

# ----------
# plot genes
# ----------

genes = wightman_hits['Gene'].unique()


for i, gene in enumerate(genes):
    print(f'Plotting gene #{i} == {gene}')
    gene_exp = expression[['id', gene]]
    gene_exp.columns = ['id', 'value']
    dkt_table = dkt_labels.merge(gene_exp, on='id')

    plot_dkt_table_brainspace(dkt_table,
                              layer='pial',
                              size=(1600, 300),
                              cmap='inferno',
                              nan_color=(0.5, 0.5, 0.5, 1),
                              filename=f'{gene}.png',
                              screenshot=True,
                              zoom=1.7)

# ----------
# Colorbar
# ----------

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

fig, ax = plt.subplots(figsize=(7, 1), dpi=300)
ax.imshow(gradient, aspect='auto', cmap='inferno')
ax.axis('off')

fig.savefig('colorbar.png', transparent=True)
