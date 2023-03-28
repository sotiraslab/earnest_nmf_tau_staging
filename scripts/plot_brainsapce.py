#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:05:59 2023

@author: tom.earnest
"""

import os

from brainstat.datasets import fetch_template_surface
from brainspace.plotting.surface_plotting import plot_hemispheres
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.io import loadmat

# find aparc.annot files
this_script = os.path.abspath(__file__)
main_dir = os.path.dirname(os.path.dirname(this_script))
resources = os.path.join(main_dir, 'resources')
LH_APARC_ANNOT = os.path.join(resources, 'lh.aparc.annot')
RH_APARC_ANNOT = os.path.join(resources, 'rh.aparc.annot')

def nmf_component_to_dkt_table(path_mat, component_index, region_names):
    mat = loadmat(path_mat)
    W = mat['Wnorm']
    component = W[:, component_index]
    dkt_table = pd.DataFrame({'region': region_names, 'value': component})

    return dkt_table

def get_vtx_data_for_hemi(dkt_table, hemi, aparc_file):

    labels, _, names = nib.freesurfer.read_annot(aparc_file)
    names = [b.decode() for b in names]

    dkt_table = dkt_table[dkt_table['region'].str.contains(f'{hemi}_')]
    region_to_value = pd.Series(data=dkt_table['value'].to_numpy(),
                                index=dkt_table['region'].str.replace('(lh|rh)_', '', regex=True))
    region_to_value = region_to_value.reindex(names)

    array_values = region_to_value.to_numpy()
    vtx_data = array_values[labels]
    vtx_data[labels == -1] = np.nan
    vtx_data[np.isnan(vtx_data)] = np.nan

    return vtx_data

def plot_dkt_table_brainspace(dkt_table, layer='inflated', **kwargs):
    ldata = get_vtx_data_for_hemi(dkt_table, 'lh', LH_APARC_ANNOT)
    rdata = get_vtx_data_for_hemi(dkt_table, 'rh', RH_APARC_ANNOT)
    data = np.concatenate([ldata, rdata])
    lsurf, rsurf = fetch_template_surface("fsaverage", join=False, layer=layer)
    fig = plot_hemispheres(
        lsurf,
        rsurf,
        data,
        **kwargs
    )
    return fig
