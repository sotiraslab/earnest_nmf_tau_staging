#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:19:42 2022

@author: tom.earnest
"""

import os
import shutil
import sys

import abagen

output_expression = 'abagen_expression_dkt.csv'
output_labels = 'dkt_labels.csv'

if os.path.exists(output_expression) and os.path.exists(output_labels):
    print("Skipping gathering of Abagen expression; "
          f"Found existing '{output_expression}' and "
          f"{output_labels}.")
    sys.exit()

dkt = abagen.datasets.fetch_desikan_killiany(surface=True)
dkt_image = dkt['image']
dkt_info = dkt['info']

# copy label info to this folder
# this raises an error, but still seems to copy correctly?
shutil.copy(dkt_info, 'dkt_labels.csv')

# get expression information
expression = abagen.allen.get_expression_data(dkt_image, verbose=2)
expression.to_csv(output_expression)
