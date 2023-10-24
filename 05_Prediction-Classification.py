#!/usr/bin/env python
# coding: utf-8



import itertools
import os
import random
import subprocess
import sys
from datetime import datetime
from functools import reduce

import openpyxl
from IPython.display import display
from tqdm import tqdm


import collections

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyreadr
import scipy as sp
from joblib import dump, load
from sklearn.datasets import make_classification
from sklearn.preprocessing import StandardScaler

from sklearn.decomposition import PCA

from sklearn.linear_model import (
    SGDClassifier, 
    SGDRegressor, 
)

from sklearn.svm import (
    SVR, SVC,
)

from sklearn.neighbors import (
    KNeighborsClassifier, 
    KNeighborsRegressor, 
)

from sklearn.ensemble import (
    ExtraTreesClassifier,
    ExtraTreesRegressor,
    HistGradientBoostingClassifier, 
    HistGradientBoostingRegressor, 
    VotingClassifier,
    VotingRegressor, 
)

from sklearn.neural_network import (
    MLPClassifier, 
    MLPRegressor
)

from sklearn.metrics import (
    RocCurveDisplay,
    accuracy_score,
    auc,
    make_scorer,
    precision_recall_curve,
    roc_auc_score,
)
from sklearn.model_selection import (
    KFold,
    cross_validate,
    train_test_split,
)


# import raw
reviewer = pyreadr.read_r(
    "XXXX"
)

DiffData = pd.read_excel("XXXX", sheet_name = None)

# data cleaning
## func
items = list(reviewer.items())
Reviewer = items[0]
Reviewer = Reviewer[1]

Reviewer_successful = Reviewer[Reviewer["successful"] == "successful"]
Reviewer_successful = Reviewer_successful.fillna(0)

TPs = ["func_TP" + str(tp) + "hr" for tp in ["0", "3", "6", "12", "24", "48", "72"]]
TPs.append("sequence_name")
func = Reviewer_successful[TPs]
dupIDs = func.duplicated()
func = func.drop_duplicates().reset_index().drop(["index"], axis=1)
TPs.remove("sequence_name")


## DiffData
sequence_name_order = func["sequence_name"].values.tolist()

Reviewer_successful = Reviewer_successful[~dupIDs]

Reviewer_successful = Reviewer_successful.set_index("sequence_name", drop=False)
Reviewer_successful = Reviewer_successful.loc[sequence_name_order]

Reviewer_successful = Reviewer_successful.reset_index(drop=True, inplace=False)

## func melt
func_melt = func.melt(
    id_vars="sequence_name", var_name="timepoint", value_name="func"
)

# "K27ac", "ATACseq", "RNAseq", "TF_RNAseq"
K27ac = Reviewer_successful[
    ["sequence_name"] + [x.replace("func_TP", "") + "_K27ac" for x in TPs]
]
ATACseq = Reviewer_successful[
    ["sequence_name"] + [x.replace("func_TP", "") + "_ATACseq" for x in TPs]
]
RNAseq = Reviewer_successful[
    ["sequence_name"] + [x.replace("func_TP", "") + "_RNAseq" for x in TPs]
]
TF_RNAseq = Reviewer_successful[
    ["sequence_name"] + [x.replace("func_TP", "") + "_TF_RNAseq" for x in TPs]
]

K27ac_melt = K27ac.melt(
    id_vars="sequence_name", var_name="timepoint", value_name="K27ac"
)
ATACseq_melt = ATACseq.melt(
    id_vars="sequence_name", var_name="timepoint", value_name="ATACseq"
)
RNAseq_melt = RNAseq.melt(
    id_vars="sequence_name", var_name="timepoint", value_name="RNAseq"
)
TF_RNAseq_melt = TF_RNAseq.melt(
    id_vars="sequence_name", var_name="timepoint", value_name="TF_RNAseq"
)

K27ac_melt.timepoint = "func_TP" + K27ac_melt.timepoint.replace("_K27ac", "", regex=True)
ATACseq_melt.timepoint = "func_TP" + ATACseq_melt.timepoint.replace("_ATACseq", "", regex=True)
RNAseq_melt.timepoint = "func_TP" + RNAseq_melt.timepoint.replace("_RNAseq", "", regex=True)
TF_RNAseq_melt.timepoint = "func_TP" + TF_RNAseq_melt.timepoint.replace("_TF_RNAseq", "", regex=True)

TPFeatures_melt = reduce(
    lambda left, right: pd.merge(left, right, on=['sequence_name', 'timepoint'], how="inner"), 
    [K27ac_melt, ATACseq_melt, RNAseq_melt, TF_RNAseq_melt]
)

func_TPFeatures_melt = TPFeatures_melt.merge(func_melt, on = ['sequence_name', "timepoint"])
func_TPFeatures_melt = func_TPFeatures_melt.set_index('sequence_name', drop = True)

# fit model
## sequence name for subsetting data in the loop
sequence_name = func[["sequence_name"]]
## TPs and Perts for loop
Perts = ["pert" + str(i) for i in list(range(1, 4))]
## correlators
Reg_correlators = ["Pearson", "Spearman", "KendalTau"]
regression_correlators = (
    (sp.stats.spearmanr, "Spearman"),
    (sp.stats.pearsonr, "Pearson"),
    (sp.stats.kendalltau, "KendalTau"),
)

# Define models with cross-validation
sgd = SGDClassifier(
    loss = 'modified_huber', 
    random_state=1, n_jobs = -1
)
svc = SVC(
    probability = True, 
    random_state=1
)
knn = KNeighborsClassifier(
    n_jobs = -1
)    
et = ExtraTreesClassifier(
    n_estimators=1500,
    random_state=1,
    n_jobs=-1,
    max_features='sqrt', 
)
gb = HistGradientBoostingClassifier(
    random_state=1,
    verbose=2,
)
mlp = MLPClassifier(
    random_state=1, 
)

models={'sgd': sgd, 
        'svc': svc, 
        'knn': knn, 
        "et": et, 
        "gb": gb, 
        "mlp": mlp}

# 10 fold cross validation
cv = KFold(n_splits=10, random_state=1, shuffle=True)

# prediction
for pertX in Perts:
    with open("log", "w") as log:
        log.write(
            "Fitting models for "
            + pertX
            + " starting "
            + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            + "\n"
        )
    #
    # get sequence ids for pertX
    ids_pertX = DiffData[pertX].sequence_name
    #
    # TP features and func
    pertX_func_TPFeatures_melt = func_TPFeatures_melt.loc[ids_pertX.tolist()]
    #
    #
    # training input --> type should be numpy array
    ## feature diff
    X_diff = DiffData[pertX].set_index(keys = "sequence_name").loc[pertX_func_TPFeatures_melt.index]# concatenating X
    ## TP features
    X_tp = pertX_func_TPFeatures_melt.drop(['timepoint', 'func'], axis = 1)
    X = pd.concat([X_diff, X_tp], axis = 1).to_numpy()
    # target values --> type should be np array
    y_func = pertX_func_TPFeatures_melt.loc[ids_pertX.tolist()]
    # sanity check
    X_tp.index.equals(y_func.index)
    # y
    y = y_func.func.to_numpy().astype('int')
    #
    for name, model in models.items():
        with open("log", "a") as log:
            log.write(
                "Fitting " 
                + name 
                + "\n"
            )
    #       
        # initiate auroc and prc
        exec("aurocs_" + pertX + "_" + name +" = []")
        #exec("auprcs_" + pertX + "_" + name +" = []")
        
        for i, (train, test) in enumerate(cv.split(X, y)):
            with open("log", "a") as log:
                log.write(
                    "Fold "
                    + str(i)
                    + " starting "
                    + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    + "\n"
                )
            #
            # fit models
            X_train = X[train]
            X_test = X[test]
            y_train = y[train]
            y_test = y[test]
            #
            model.fit(X_train, y_train) 
            #
            # predicting
            y_pred = model.predict_proba(X_test)
            #
            ## calc AUROC
            auroc = roc_auc_score(y_test, y_pred, multi_class = "ovr")
            #
            ## calc AUPRC
            #precision, recall, thresholds = precision_recall_curve(y_test, y_pred[:, 1])
            #auprc = auc(recall, precision)
            #
            exec("aurocs_" + pertX + "_" + name +".append(auroc)")
            #exec("auprcs_" + pertX + "_" + name +".append(auprc)")


for name, model in models.items():
    exec("aurocs_" + name 
         + " = pd.DataFrame([aurocs_pert1_" + name 
         + ", aurocs_pert2_" + name 
         + ", aurocs_pert3_" + name 
         + "], index=[f'pert {i}' for i in range(1, 4)], columns=[f'Fold {i}' for i in range(1, 11)],).transpose()")
    exec("aurocs_" + name + ".to_csv('aurocs_' + name + '.csv')")
