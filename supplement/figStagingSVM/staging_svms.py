#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 14:49:32 2024

@author: earnestt1234
"""

#%% Imports

import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import LinearSVC
from scipy.stats import t

#%% Paths

DATA_ADNI = '../../derivatives/adni/data_with_staging.csv'
DATA_OASIS = '../../derivatives/oasis3/data_with_staging.csv'

#%% Read

adni = pd.read_csv(DATA_ADNI)
oasis = pd.read_csv(DATA_OASIS)

#%% Create data for experiements

# omit missing
adni.dropna(subset=['CDRBinned', 'Age', 'Gender', 'Centiloid', 'HasE4'], inplace=True)
oasis.dropna(subset=['CDR', 'Age', 'Gender', 'Centiloid', 'HasE4'], inplace=True)

# make equivalent CDR columns
adni['CDRTarget'] = adni['CDRBinned'].copy()
adni = adni.loc[adni['CDRTarget'] != '']

oasis['CDRTarget'] = oasis['CDR'].astype(str)
oasis.loc[oasis['CDRTarget'] == '1.0', 'CDRTarget'] = '1.0+'


# separate
adni_ads = adni[adni['Group'] == 'TrainingBaseline'].copy()
oasis_ads = oasis[oasis['Group'] == 'TrainingSet'].copy()

#%% Helper functions

def make_feature_matrix(dataset, feature_cols):
    tmp = []
    for f in feature_cols:
        data = dataset[f]
        if is_numeric_dtype(data):
            tmp.append(data.to_numpy()[:, np.newaxis])
        else:
            tmp.append(pd.get_dummies(data, drop_first=True))

    X = np.hstack(tmp)
    return X

def svm_experiment(adni_data, oasis_data, feature_cols,
                   target_col='CDRTarget', stratify_col='CDRTarget',
                   n_repeats=10, n_folds=10,
                   metric=accuracy_score):

    internal = adni_data
    external = oasis_data

    # create features
    X_internal = make_feature_matrix(internal, feature_cols)
    X_external = make_feature_matrix(external, feature_cols)

    # create target
    target_internal = pd.factorize(internal[target_col])[0]
    target_external = pd.factorize(external[target_col])[0]

    # main loop: model training
    internal_scores = []
    models = []
    sizes_train = []
    sizes_test = []

    for r in range(n_repeats):

        cv = StratifiedKFold(n_splits=n_folds, random_state=r, shuffle=True)

        for i, (train_index, test_index) in enumerate(cv.split(internal, internal[stratify_col])):

            # ptc model
            train_X = X_internal[train_index, :]
            test_X = X_internal[test_index, :]

            train_target = target_internal[train_index]
            test_target = target_internal[test_index]

            model = LinearSVC(C=1, class_weight='balanced', dual=False)
            model.fit(train_X, train_target)
            preds = model.predict(test_X)
            score = metric(test_target, preds)
            internal_scores.append(score)
            models.append(model)
            sizes_train.append(len(train_X))
            sizes_test.append(len(test_X))

    # external prediction
    external_scores = [metric(target_external, m.predict(X_external)) for m in models]

    return {'adni_scores': internal_scores,
            'oasis_scores': external_scores,
            'models': models,
            'sizes_train': sizes_train,
            'sizes_test': sizes_test}

def nadeau_bengio_test(a, b, n_train, n_test, alpha=0.05, side='both'):
    """
    Implementation of the Nadeau Bengio correction
    for the t-test [1].  This is recommended for accounting
    for the dependence of data points when comparing cross-validated
    model performance.

    Formula is based on the equation outlined by
    Bouckaert & Frank [2]; see section 3.2.

    Implementation also follows this Gist [3], but that
    has some errors that are fixed here.

    [1] https://proceedings.neurips.cc/paper/1999/hash/7d12b66d3df6af8d429c1a357d8b9e1a-Abstract.html
    [2] https://link.springer.com/chapter/10.1007/978-3-540-24775-3_3
    [3] https://gist.github.com/jensdebruijn/13e8eeda85eb8644ac2a4ac4c3b8e732

    Parameters
    ----------
    a : array
        Model A performance metrics.
    b : array
        Model B performance metrics.
    n_train : int
        Number of observations used in a single training fold.
    n_test : int
        Number of observations used in a single testing fold.

    Returns
    -------
    None.

    """

    # check arguments
    if len(a) != len(b):
        raise ValueError("`a` and `b` inputs must be arrays of same length.")

    if side not in ['left', 'right', 'both']:
        raise ValueError("`side` must be 'left', 'right', or 'both'")

    # set variables for equation
    x = np.array(a) - np.array(b)
    var = np.var(x, ddof=1)
    n = len(x)
    n2 = n_test
    n1 = n_train

    # calculate statistic
    numerator = np.mean(x)
    parens = ((1/n) + (n2/n1))
    denominator = np.sqrt(parens * var)
    tstat = numerator / denominator

    # calculate p-value
    dof = n - 1

    if side == 'left':
        p = t.cdf(tstat, dof)
    elif side == 'right':
        p = 1 - t.cdf(tstat, dof)
    elif side == 'both':
        p = 2 * (1 - t.cdf(abs(tstat), dof))

    return {'t': tstat, 'p': p, 'dof': dof}

def compare_experiments(a, b):
    n_train = int(np.mean(a['sizes_train']))
    n_test = int(np.mean(a['sizes_test']))

    # adni
    ttest = nadeau_bengio_test(a['adni_scores'], b['adni_scores'],
                               n_train=n_train, n_test=n_test,
                               alpha=0.05, side='both')
    adni_stats = {'dataset': 'adni',
                  'mean_a': np.mean(a['adni_scores']),
                  'std_a': np.std(a['adni_scores']),
                  'mean_b': np.mean(b['adni_scores']),
                  'std_b': np.std(b['adni_scores']),
                  't': ttest['t'],
                  'p': ttest['p'],
                  'dof': ttest['dof']}

    # oasis
    ttest = nadeau_bengio_test(a['oasis_scores'], b['oasis_scores'],
                                     n_train=n_train, n_test=n_test,
                                     alpha=0.05, side='both')
    oasis_stats = {'dataset': 'oasis',
                   'mean_a': np.mean(a['oasis_scores']),
                   'std_a': np.std(a['oasis_scores']),
                   'mean_b': np.mean(b['oasis_scores']),
                   'std_b': np.std(b['oasis_scores']),
                   't': ttest['t'],
                   'p': ttest['p'],
                   'dof': ttest['dof']}

    return adni_stats, oasis_stats

#%% Run models

# staging only models
svm_ptc = svm_experiment(adni_ads, oasis_ads, feature_cols=['PTCStage'])
svm_braak= svm_experiment(adni_ads, oasis_ads, feature_cols=['BraakStage'])

# other covariates included
svm_ptc_cov = svm_experiment(adni_ads, oasis_ads, feature_cols=['PTCStage', 'Age', 'Gender', 'HasE4', 'Centiloid'])
svm_braak_cov = svm_experiment(adni_ads, oasis_ads, feature_cols=['BraakStage', 'Age', 'Gender', 'HasE4', 'Centiloid'])

#%% Collect results

pairs = [
    [svm_ptc, svm_braak],
    [svm_ptc_cov, svm_braak_cov],
    ]

covariates = [False, True]

rows = []
for i, pair in enumerate(pairs):
    adni_results, oasis_results = compare_experiments(*pair)
    cov = 'Age, Sex, APOE E4+, Centiloid' if covariates[i] else '-'
    rows.append({'Dataset': 'ADNI',
                 'Covariates': cov,
                 'PTC': f"{round(adni_results['mean_a'], 2)} ({(round(adni_results['std_a'], 3))})",
                 'Braak': f"{round(adni_results['mean_b'], 2)} ({(round(adni_results['std_b'], 3))})",
                 'p-value': round(adni_results['p'], 3),
                 't-value': round(adni_results['t'], 3)})

    rows.append({'Dataset': 'OASIS',
                 'Covariates': cov,
                 'PTC': f"{round(oasis_results['mean_a'], 2)} ({(round(oasis_results['std_a'], 3))})",
                 'Braak': f"{round(oasis_results['mean_b'], 2)} ({(round(oasis_results['std_b'], 3))})",
                 'p-value': round(oasis_results['p'], 3),
                 't-value': round(oasis_results['t'], 3)})


output = pd.DataFrame(rows).sort_values(['Dataset'])
output.to_csv('svm_results_formatted.csv', index=False)
