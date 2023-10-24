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
## Log2FC
items = list(reviewer.items())
Reviewer = items[0]
Reviewer = Reviewer[1]

Reviewer_successful = Reviewer[Reviewer["successful"] == "successful"]
Reviewer_successful = Reviewer_successful.fillna(0)

TPs = ["Log2FC_TP" + str(tp) + "hr" for tp in ["0", "3", "6", "12", "24", "48", "72"]]
TPs.append("sequence_name")
Log2FC = Reviewer_successful[TPs]
dupIDs = Log2FC.duplicated()
Log2FC = Log2FC.drop_duplicates().reset_index().drop(["index"], axis=1)
TPs.remove("sequence_name")


## DiffData
sequence_name_order = Log2FC["sequence_name"].values.tolist()

Reviewer_successful = Reviewer_successful[~dupIDs]

Reviewer_successful = Reviewer_successful.set_index("sequence_name", drop=False)
Reviewer_successful = Reviewer_successful.loc[sequence_name_order]

Reviewer_successful = Reviewer_successful.reset_index(drop=True, inplace=False)

## Log2FC melt
Log2FC_melt = Log2FC.melt(
    id_vars="sequence_name", var_name="timepoint", value_name="Log2FC"
)

# "K27ac", "ATACseq", "RNAseq", "TF_RNAseq"
K27ac = Reviewer_successful[
    ["sequence_name"] + [x.replace("Log2FC_TP", "") + "_K27ac" for x in TPs]
]
ATACseq = Reviewer_successful[
    ["sequence_name"] + [x.replace("Log2FC_TP", "") + "_ATACseq" for x in TPs]
]
RNAseq = Reviewer_successful[
    ["sequence_name"] + [x.replace("Log2FC_TP", "") + "_RNAseq" for x in TPs]
]
TF_RNAseq = Reviewer_successful[
    ["sequence_name"] + [x.replace("Log2FC_TP", "") + "_TF_RNAseq" for x in TPs]
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

K27ac_melt.timepoint = "Log2FC_TP" + K27ac_melt.timepoint.replace("_K27ac", "", regex=True)
ATACseq_melt.timepoint = "Log2FC_TP" + ATACseq_melt.timepoint.replace("_ATACseq", "", regex=True)
RNAseq_melt.timepoint = "Log2FC_TP" + RNAseq_melt.timepoint.replace("_RNAseq", "", regex=True)
TF_RNAseq_melt.timepoint = "Log2FC_TP" + TF_RNAseq_melt.timepoint.replace("_TF_RNAseq", "", regex=True)

TPFeatures_melt = reduce(
    lambda left, right: pd.merge(left, right, on=['sequence_name', 'timepoint'], how="inner"), 
    [K27ac_melt, ATACseq_melt, RNAseq_melt, TF_RNAseq_melt]
)

Log2FC_TPFeatures_melt = TPFeatures_melt.merge(Log2FC_melt, on = ['sequence_name', "timepoint"])
Log2FC_TPFeatures_melt = Log2FC_TPFeatures_melt.set_index('sequence_name', drop = True)

# fit model
## sequence name for subsetting data in the loop
sequence_name = Log2FC[["sequence_name"]]
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
sgd = SGDRegressor(
    loss = 'huber', 
    random_state=1, 
)
svr = SVR()
knn = KNeighborsRegressor(
    n_jobs = -1
)    
et = ExtraTreesRegressor(
    n_estimators=1500,
    random_state=1,
    n_jobs=-1,
    max_features='sqrt', 
)
gb = HistGradientBoostingRegressor(
    random_state=1,
    verbose=2,
)
mlp = MLPRegressor(
    random_state=1, 
)

# 10 fold cross validation
cv = KFold(n_splits=10, random_state=1, shuffle=True)

# initiate empty data frames for correlation coeficients and p values
CorCoef = pd.DataFrame(np.nan, index=list(range(0, 10)), columns=Reg_correlators)
CorCoef.columns.name = "correlator"
CorCoef.index.name = "Fold"
CorP = CorCoef.copy()

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
    # TP features and Log2FC
    pertX_Log2FC_TPFeatures_melt = Log2FC_TPFeatures_melt.loc[ids_pertX.tolist()]
    #
    #
    # training input --> type should be numpy array
    ## feature diff
    X_diff = DiffData[pertX].set_index(keys = "sequence_name").loc[pertX_Log2FC_TPFeatures_melt.index]# concatenating X
    ## TP features
    X_tp = pertX_Log2FC_TPFeatures_melt.drop(['timepoint', 'Log2FC'], axis = 1)
    X = pd.concat([X_diff, X_tp], axis = 1).to_numpy()
    # target values --> type should be np array
    y_Log2FC = pertX_Log2FC_TPFeatures_melt.loc[ids_pertX.tolist()]
    # sanity check
    X_tp.index.equals(y_Log2FC.index)
    # y
    y = y_Log2FC.Log2FC.to_numpy()#.astype('int')
    #
    models={'sgd': sgd, 
        'svr': svr, 
        'knn': knn, 
        "et": et, 
        "gb": gb, 
        "mlp": mlp}
    for name, model in models.items():
        with open("log", "a") as log:
            log.write(
                "Fitting " 
                + name 
                + "\n"
            )
        #
        # initiate auroc and prc
        exec("CorCoef_" + pertX + "_" + name +" = CorCoef.copy()")
        exec("CorP_" + pertX + "_" + name +" = CorP.copy()")
        #  
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
            y_pred = model.predict(X_test)
            #
            for correlator, label in regression_correlators:
                exec("CorCoef_" + pertX + "_" + name + ".loc[i][label] = correlator(y_pred, y_test)[0]")
                exec("CorP_" + pertX + "_" + name + ".loc[i][label] = correlator(y_pred, y_test)[1]")


for name, model in models.items(): 
    exec("writerCoef = pd.ExcelWriter('CorCoef_' + name + '.xlsx', mode = 'w')")
    for pertX in Perts: 
        exec('CorCoef_' + pertX + '_' + name + '.to_excel(writerCoef, sheet_name = "'+ pertX + '", index=False, engine="xlsxwriter")')
        writerCoef.save()

for name, model in models.items(): 
    exec("writerP = pd.ExcelWriter('CorP_' + name + '.xlsx', mode = 'w')")
    for pertX in Perts: 
        exec('CorP_' + pertX + '_' + name + '.to_excel(writerP, sheet_name = "'+ pertX + '", index=False, engine="xlsxwriter")')
        writerP.save()
