#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 12:25:59 2023

@author: earnestt1234
"""

import os
import re

import matplotlib.pyplot as plt
import natsort
import numpy as np
from scipy.io import loadmat

def recon_error_plots(input_matrix_csv, nmf_results_folder,
                      output_folder=None, scale=False,
                      xlab='Number of Components'):

    # read data
    X = np.loadtxt(input_matrix_csv, delimiter=',')

    errors = []
    ranks = []
    for mat in natsort.natsorted(os.listdir(nmf_results_folder)):
        if not mat.lower().endswith('mat'): continue;
        print(f'{mat}...')
        fullfile = os.path.join(nmf_results_folder, mat)
        d = loadmat(fullfile)
        W = d['W']
        H = d['H']
        fro = np.linalg.norm(X - np. matmul(W, H), ord='fro')
        errors.append(fro)
        ranks.append(int(re.search('\d+', mat).group()))

    errors = np.array(errors)
    if scale:
        errors = (errors - errors.min()) / (errors.max() - errors.min())
    # ranks = np.array(range(2, 21))

    # plots
    plt.rcParams.update({'font.family':'arial',
                         'font.size': 24})

    # recon error
    xticks = [i for i in ranks if i % 2 == 0]
    ylab = 'Reconstruction error'

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(ranks, errors, zorder=3, color='red', lw=3)

    ax.set_xticks(xticks)
    ax.grid(zorder=1, alpha=.3)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)

    if output_folder is not None:
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, 'recon_error.png'), dpi=300)

    # gradient recon error
    ylab = 'Gradient of reconstruction error'

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(ranks, np.diff(errors, append=np.nan), zorder=3, color='red', lw=3)

    ax.set_xticks(xticks)
    ax.grid(zorder=1, alpha=.3)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)

    if output_folder is not None:
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, 'recon_error_gradient.png'), dpi=300)

def reproducibility_plots(reproducibility_stats_mat, output_folder=None,
                          xlab='Number of Components'):

    d = loadmat(reproducibility_stats_mat)
    ari = d['repeat_ARI']
    meanInner = d['repeat_meanInner']
    medianInner = d['repeat_medianInner']
    ranks = d['sortedBasisNum'].flatten()
    xticks = [i for i in ranks if i % 2 == 0]

    # plots
    plt.rcParams.update({'font.family':'arial',
                         'font.size': 24})

    # ari mean
    x = ranks

    color = 'dodgerblue'
    ylab = 'Adjusted Rand Index'

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x, ari.mean(axis=0), zorder=3, lw=3, color=color)

    ax.set_xticks(xticks)
    ax.grid(zorder=1, alpha=.3)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)

    if output_folder is not None:
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, 'reproducibility_ari.png'), dpi=300)

    # mean inner product
    x = ranks
    color = 'dodgerblue'
    ylab = 'Mean Inner Product'

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(x, meanInner.mean(axis=0), zorder=3, lw=3, color=color)
    ax.set_xticks(xticks)
    ax.grid(zorder=1, alpha=.3)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)

    if output_folder is not None:
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, 'reproduciblity_mean_inner.png'), dpi=300)
