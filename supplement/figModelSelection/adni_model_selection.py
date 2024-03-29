#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:22:24 2023

@author: earnestt1234
"""

import sys

# paths
PATH_SCRIPTS = '../../scripts/'
PATH_INPUT_MAT = '../../derivatives/adni/nmf_matrix.csv'
PATH_RESULTS = '../../nmf/adni/results/mat'
PATH_REPRODUCIBILITY = '../../nmf/adni/results/repeat_reproducibility_mat/reproducibilityResults.mat'

# load NMF model selection plot code
sys.path.append(PATH_SCRIPTS)
from nmf_model_selection_plots import recon_error_plots, reproducibility_plots

# recon error
recon_error_plots(PATH_INPUT_MAT, PATH_RESULTS, output_folder='.',
                  xlab='Number of PTCs')


# reproducibility_plots
reproducibility_plots(PATH_REPRODUCIBILITY, output_folder='.',
                      xlab='Number of PTCs')
