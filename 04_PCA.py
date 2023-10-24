#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sklearn.decomposition import PCA

import pyreadr
import collections
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import seaborn as sns

from sklearn.preprocessing import StandardScaler


# # import raw data


Data = "XXXX"

object_list = pyreadr.list_objects(
    "XXXX"
)

reviewer = pyreadr.read_r(
    "XXXX"
)

DiffData = pd.read_csv(Data + "diff.csv")

DiffData = pd.read_csv(Data + "diff.csv")


# # data cleaning
# 
# ## Reviewer and Log2FC

items = list(reviewer.items())
Reviewer = items[0]
Reviewer = Reviewer[1]

Reviewer_successful = Reviewer[Reviewer["successful"] == "successful"]
Reviewer_successful = Reviewer_successful.fillna(0)


Reviewer_successful.drop_duplicates().shape

TPs = ["Log2FC_TP" + str(tp) + "hr" for tp in ["0", "3", "6", "12", "24", "48", "72"]]
TPs.append("sequence_name")
Log2FC = Reviewer_successful[TPs]
dupIDs = Log2FC.duplicated()
Log2FC = Log2FC.drop_duplicates().reset_index().drop(["index"], axis=1)
TPs.remove("sequence_name")


# ## Reviewer successful

sequence_name_order = Log2FC["sequence_name"].values.tolist()
Reviewer_successful = Reviewer_successful[~dupIDs]
Reviewer_successful = Reviewer_successful.set_index("sequence_name", drop=False)
Reviewer_successful = Reviewer_successful.loc[sequence_name_order]

Reviewer_successful = Reviewer_successful.reset_index(drop=True, inplace=False)


# #  Sanity check

## sanity check
Log2FC[["sequence_name"]].equals(DiffData[["sequence_name"]])
DiffData[["sequence_name"]].equals(pd.DataFrame(Reviewer_successful.index))
# ### prepare diff data for fitting models
DiffData_features = DiffData.drop(["sequence_name"], axis=1)
DiffData_features = DiffData.set_index(keys = "sequence_name")


# # PCA

# TPs and Perts for loop
Perts = ["pert" + str(i) for i in list(range(1, 4))]


writerPCACoord = pd.ExcelWriter("PCACoord.xlsx")
writerPCALoading = pd.ExcelWriter("PCALoading.xlsx")


for pertX in Perts:
    with open("PCA.log", "w") as log:
        log.write(
            "PCA models for "
            + pertX
            + " starting "
            + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            + "\n"
        )
    # get Diff
    ## get sequence ids for pertX
    ids_pertX = DiffData_features.index[DiffData_features.index.str.contains(pertX)]
    ## feature diff
    X_diff = DiffData_features.loc[ids_pertX]# concatenating X
    #
    # scale
    sc = StandardScaler()
    sc.fit(X_diff)
    X_diff_std = sc.transform(X_diff)
    #
    # PCA
    pca = PCA(n_components = 0.99)
    pca_df = pca.fit_transform(X_diff_std)
    pca_coord = pd.DataFrame(pca_df)
    pca_coord = pca_coord.set_index(keys = ids_pertX)
    #
    # plot PCA
    # Determine explained variance using explained_variance_ration_ attribute
    exp_var_pca = pca.explained_variance_ratio_
    #
    # Cumulative sum of eigenvalues; This will be used to create step plot
    # for visualizing the variance explained by each principal component.
    #
    cum_sum_eigenvalues = np.cumsum(exp_var_pca)
    #
    # Create the visualization plot
    #
    plt.bar(range(0,len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
    plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal component index')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.title(pertX)
    plt.show()
    #
    loadings = pca.components_
    num_pc = pca.n_features_
    pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
    loadings_X_diff = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
    loadings_X_diff['variable'] = X_diff.columns.values
    loadings_X_diff = loadings_X_diff.set_index('variable')
    #
    # write out 
    pca_coord.to_excel(writerPCACoord, pertX, index=True, engine="xlsxwriter")
    writerPCACoord.save()
    loadings_X_diff.to_excel(writerPCALoading, pertX, index=True, engine="xlsxwriter")
    writerPCALoading.save()    


pertX = "pert1"
with open("PCA.log", "w") as log:
    log.write(
        "PCA models for "
        + pertX
        + " starting "
        + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        + "\n"
    )
# get Diff
## get sequence ids for pertX
ids_pertX = DiffData_features.index[DiffData_features.index.str.contains(pertX)]
## feature diff
X_diff = DiffData_features.loc[ids_pertX]# concatenating X
#
# scale
sc = StandardScaler()
sc.fit(X_diff)
X_diff_std = sc.transform(X_diff)
#
# PCA
pca = PCA(n_components = 0.9999)
pca_df = pca.fit_transform(X_diff_std)
pca_coord = pd.DataFrame(pca_df)
pca_coord = pca_coord.set_index(keys = ids_pertX)
#
# plot PCA
# Determine explained variance using explained_variance_ration_ attribute
exp_var_pca = pca.explained_variance_ratio_
#
# Cumulative sum of eigenvalues; This will be used to create step plot
# for visualizing the variance explained by each principal component.
#
cum_sum_eigenvalues = np.cumsum(exp_var_pca)
#
# Create the visualization plot
#
plt.bar(range(0,len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')
plt.ylabel('Explained variance ratio')
plt.xlabel('Principal component index')
plt.legend(loc='best')
plt.tight_layout()
plt.title(pertX)
plt.show()
#
loadings = pca.components_
num_pc = pca.n_features_
pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
loadings_X_diff = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
loadings_X_diff['variable'] = X_diff.columns.values
loadings_X_diff = loadings_X_diff.set_index('variable')
#
# write out 
pca_coord.to_excel(writerPCACoord, pertX, index=True, engine="xlsxwriter")
writerPCACoord.save()
#loadings_X_diff.to_excel(writerPCALoading, pertX, index=True, engine="xlsxwriter")
#writerPCALoading.save()
loadings_X_diff.to_csv("./PCALoading_" + pertX + ".csv")


pertX = "pert2"
with open("PCA.log", "w") as log:
    log.write(
        "PCA models for "
        + pertX
        + " starting "
        + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        + "\n"
    )
# get Diff
## get sequence ids for pertX
ids_pertX = DiffData_features.index[DiffData_features.index.str.contains(pertX)]
## feature diff
X_diff = DiffData_features.loc[ids_pertX]# concatenating X
#
# scale
sc = StandardScaler()
sc.fit(X_diff)
X_diff_std = sc.transform(X_diff)
#
# PCA
pca = PCA(n_components = 0.9999)
pca_df = pca.fit_transform(X_diff_std)
pca_coord = pd.DataFrame(pca_df)
pca_coord = pca_coord.set_index(keys = ids_pertX)
#
# plot PCA
# Determine explained variance using explained_variance_ration_ attribute
exp_var_pca = pca.explained_variance_ratio_
#
# Cumulative sum of eigenvalues; This will be used to create step plot
# for visualizing the variance explained by each principal component.
#
cum_sum_eigenvalues = np.cumsum(exp_var_pca)
#
# Create the visualization plot
#
plt.bar(range(0,len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')
plt.ylabel('Explained variance ratio')
plt.xlabel('Principal component index')
plt.legend(loc='best')
plt.tight_layout()
plt.title(pertX)
plt.show()
#
loadings = pca.components_
num_pc = pca.n_features_
pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
loadings_X_diff = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
loadings_X_diff['variable'] = X_diff.columns.values
loadings_X_diff = loadings_X_diff.set_index('variable')
#
# write out 
pca_coord.to_excel(writerPCACoord, pertX, index=True, engine="xlsxwriter")
writerPCACoord.save()
#loadings_X_diff.to_excel(writerPCALoading, pertX, index=True, engine="xlsxwriter")
#writerPCALoading.save()
loadings_X_diff.to_csv("./PCALoading_" + pertX + ".csv")


pertX = "pert3"
with open("PCA.log", "w") as log:
    log.write(
        "PCA models for "
        + pertX
        + " starting "
        + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        + "\n"
    )
# get Diff
## get sequence ids for pertX
ids_pertX = DiffData_features.index[DiffData_features.index.str.contains(pertX)]
## feature diff
X_diff = DiffData_features.loc[ids_pertX]# concatenating X
#
# scale
sc = StandardScaler()
sc.fit(X_diff)
X_diff_std = sc.transform(X_diff)
#
# PCA
pca = PCA(n_components = 0.9999)
pca_df = pca.fit_transform(X_diff_std)
pca_coord = pd.DataFrame(pca_df)
pca_coord = pca_coord.set_index(keys = ids_pertX)
#
# plot PCA
# Determine explained variance using explained_variance_ration_ attribute
exp_var_pca = pca.explained_variance_ratio_
#
# Cumulative sum of eigenvalues; This will be used to create step plot
# for visualizing the variance explained by each principal component.
#
cum_sum_eigenvalues = np.cumsum(exp_var_pca)
#
# Create the visualization plot
#
plt.bar(range(0,len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')
plt.ylabel('Explained variance ratio')
plt.xlabel('Principal component index')
plt.legend(loc='best')
plt.tight_layout()
plt.title(pertX)
plt.show()
#
loadings = pca.components_
num_pc = pca.n_features_
pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
loadings_X_diff = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
loadings_X_diff['variable'] = X_diff.columns.values
loadings_X_diff = loadings_X_diff.set_index('variable')
#
# write out 
pca_coord.to_excel(writerPCACoord, pertX, index=True, engine="xlsxwriter")
writerPCACoord.save()
#loadings_X_diff.to_excel(writerPCALoading, pertX, index=True, engine="xlsxwriter")
#writerPCALoading.save()
loadings_X_diff.to_csv("./PCALoading_" + pertX + ".csv")
