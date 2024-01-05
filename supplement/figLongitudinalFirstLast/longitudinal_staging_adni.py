#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 10:06:26 2022

@author: earnestt1234
"""

import pandas as pd

df = pd.read_csv('../../derivatives/adni/data_with_staging.csv')

subject_col = 'RID'
position_col = 'VisitNumber'
label_col = 'PTCStage'

df = df[[subject_col, position_col, label_col]]

subjects_vc = df[subject_col].value_counts()
subjects_long = subjects_vc.index[subjects_vc > 1]
df = df.loc[df[subject_col].isin(subjects_long)]

labels = sorted(df[label_col].unique(), key = lambda x: (x.isnumeric(), x))
positions = sorted(df[position_col].unique())

g = df.groupby(subject_col)

# CHANGE HERE!
# Group to select first and last of each subject
# rather than all visits
flow = pd.DataFrame({'subject': g[subject_col].first(),
                     'pos_from': g[position_col].first(),
                     'pos_to': g[position_col].last(),
                     'label_from': g[label_col].first(),
                     'label_to': g[label_col].last()})

flow.to_csv('adni_stage_flow.csv', index=False)

#%% report some stats

def label(x):

    if (x['label_to'] == 'NS') or (x['label_from'] == 'NS'):
        return 'NS'
    elif int(x['label_to']) == int(x['label_from']):
        return 'stable'
    elif int(x['label_to']) > int(x['label_from']):
        return 'increasing'
    elif int(x['label_to']) < int(x['label_from']):
        return 'decreasing'
    else:
        return None

v = flow.apply(label, axis=1)
print(v.value_counts() / len(v))
