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

# remove NS
# df = df[df[label_col] != 'NS']

subjects_vc = df[subject_col].value_counts()
subjects_long = subjects_vc.index[subjects_vc > 1]
df = df.loc[df[subject_col].isin(subjects_long)]


labels = sorted(df[label_col].unique(), key = lambda x: (x.isnumeric(), x))
positions = sorted(df[position_col].unique())

circle_sizes = df.groupby([position_col, label_col]).size().rename('Size')
max_circle = circle_sizes.max()
min_circle = circle_sizes.min()
circle_sizes = (circle_sizes.
                reindex(
                    pd.MultiIndex.from_product(
                        [positions, labels],
                        names=[position_col, label_col])
                    ).
                fillna(0).
                reset_index()
                )


g = df.groupby(subject_col)
flow = pd.DataFrame({'subject': df[subject_col],
                     'pos_from': df[position_col],
                     'pos_to': g[position_col].shift(-1),
                     'label_from': df[label_col],
                     'label_to': g[label_col].shift(-1)}).dropna()

flow.to_csv('adni_stage_flow.csv', index=False)

#%% report som stats

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



# a = flow['label_from']

# b = flow['label_to']
# n = len(flow)

# print()
# print('Stable: {}'.format((a == b).sum() / n))
# print('Any NS: {}'.format((a )))
