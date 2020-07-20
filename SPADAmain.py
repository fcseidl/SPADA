#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 16:07:51 2020

@author: fcseidl

This module contains various tests.
"""

import numpy as np
from sklearn.decomposition import PCA, FactorAnalysis
import matplotlib.pyplot as plt
import csv
from sklearn.preprocessing import normalize
from scipy import stats

from ZIFA import ZIFA

import preprocessing
import simulations as sims
import HypothesisTesting as ht
import SPADAutil as util


# CH relatedness testing on simulated scRNA-seq datasets
if 1:
    n_genes = 273
    n_types = 3
    lam = 0.1
    Lbig = 213
    Lsmall = 90
    n_components = 10   # TODO: choose this dynamically?
    
    print("Generating joint scRNA-seq datasets Y1 and Y2, and unrelated dataset Y3...")
    alpha1 = sims.randomalpha(n_types)
    A1 = sims.randomA(n_genes, n_types)
    A2 = sims.randomA(n_genes, n_types)
    
    _, Y1 = sims.simulateJointData(L=Lbig, A=A1, alpha=alpha1, lam=lam)
    _, Y2 = sims.simulateJointData(L=Lsmall, A=A1, alpha=alpha1, lam=lam)
    _, Y3 = sims.simulateJointData(L=Lsmall, A=A2, lam=lam) # unrelated
    
    print("\nPreprocessing Y1 and Y2 for ZIFA...")
    Y12 = np.concatenate((Y1, Y2), axis=1)  # joined data matrices
    Y12 = preprocessing.ZIFApreprocessing(Y12)
    
    print("\nPreprocessing Y1 and Y3 for ZIFA...")
    Y13 = np.concatenate((Y1, Y3), axis=1)
    Y13 = preprocessing.ZIFApreprocessing(Y13)
    
    print("\nPerforming ZIFA on data Y1 and Y2...")
    Y12, _ = ZIFA.fitModel(Y12.T, n_components)
    Y12 = Y12.T             # dim-reduced joined data matrix
    Y1 = Y12[:, :Lbig]      # columns corresponding to Y1
    Y2 = Y12[:, -Lsmall:]   # columns corresponding to Y2
    
    # TODO: avoid using n_types, hidden info
    print("\nCluster heterogeneity on Y1 and Y2 (expect high heterogeneity):")
    ht.clusterHeterogeneity(Y1, Y2)#, n_clusters=n_types)
    
    print("\nPerforming ZIFA on data Y1 and Y3...")
    Y13, _ = ZIFA.fitModel(Y13.T, n_components)
    Y13 = Y13.T
    Y1 = Y13[:, :Lbig]
    Y3 = Y13[:, -Lsmall:]
    
    print("\nCH on Y1 and Y3 (expect low heterogeneity):")
    ht.clusterHeterogeneity(Y1, Y3)#, n_clusters=2*n_types)

# how much variance is explained by principal components?
if 0:
    '''
    print("Dataset: pancreatic islets single cell")
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_sc.csv"
    '''
    print("Dataset: 3cl single cell")
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    
    print("Reading data matrix Y...")
    Y = preprocessing.csvToMatrix(scfile)
    
    print("Performing PCA...")
    pca = PCA(n_components=0.9)   # Explain at least 90% of variance
    pca.fit(Y.T)
    
    print("Inferred dimensionality:", pca.components_.shape[0])
    print("Explained percentages of variance:", pca.explained_variance_ratio_)

# hypothesis testing with 2 real datasets
if 0:
    np.random.seed(0)
    
    print("Bulk: pancreatic islets")
    print("Single-cell: pancreatic islets")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()
    
    print("Bulk: 3cl mixture")
    print("Single-cell: 3cl mixture")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()
    
    print("Bulk: 3cl mixture")
    print("Single-cell: pancreatic islets")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_islets_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_islets_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()
    
    print("Bulk: pancreatic islets")
    print("Single-cell: 3cl mixture")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_3cl_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_3cl_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()

# hypothesis testing with simulated data
if 0:
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X1, Y1 = sims.simulateJointData()
    X2, Y2 = sims.simulateJointData()
    
    '''
    print("Log-tranforming...")
    X1 = np.log(X1 + 1)
    X2 = np.log(X2 + 1)
    Y1 = np.log(Y1 + 1)
    '''
    
    '''
    print("Monotonically tranforming...")
    c = 10
    X1 = sims.g_star(X1, c)
    X2 = sims.g_star(X2, c)
    Y1 = sims.g_star(Y1, c)
    '''
    
    print('Are X1 and Y1 joint? Expect 0, receive', 
          ht.residualsToCone(X1, Y1))
    print('Are X2 and Y1 joint? Expect > 0, receive', 
          ht.residualsToCone(X2, Y1))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          ht.pvalue(X1, Y1))
    print("Probability of X2 and Y1 under null hypothesis <=",
          ht.pvalue(X2, Y1))

# hypothesis testing on dimensionality-reduced simulated data with PCA
if 0:
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X1, Y1 = sims.simulateJointData()
    X2, Y2 = sims.simulateJointData()
    
    print("Performing PCA on single-cell data...")
    D = 0.999999
    pca = PCA(n_components=D)
    Y1_hat = pca.fit_transform(Y1.T).T
    
    print("Transforming bulk data into feature space...")
    X1_hat = pca.transform(X1.T).T
    X2_hat = pca.transform(X2.T).T
    
    print('Are X1 and Y1 joint? Expect 0, receive', 
          ht.residualsToCone(X1_hat, Y1_hat))
    print('Are X2 and Y1 joint? Expect > 0, receive', 
          ht.residualsToCone(X2_hat, Y1_hat))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          ht.pvalue(X1_hat, Y1_hat))
    print("Probability of X2 and Y1 under null hypothesis <=",
          ht.pvalue(X2_hat, Y1_hat))
    
# hypothesis testing on dimensionality-reduced simulated data with FA
if 0:
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X1, Y1 = sims.simulateJointData()
    X2, Y2 = sims.simulateJointData()
    
    print("Performing Factor Analysis on single-cell data...")
    D = 25
    fa = FactorAnalysis(n_components=D)
    Y1_hat = fa.fit_transform(Y1.T).T
    
    print("Transforming bulk data into feature space...")
    X1_hat = fa.transform(X1.T).T
    X2_hat = fa.transform(X2.T).T
    
    print('Are X1 and Y1 joint? Expect 0, receive', 
          ht.residualsToCone(X1_hat, Y1_hat))
    print('Are X2 and Y1 joint? Expect > 0, receive', 
          ht.residualsToCone(X2_hat, Y1_hat))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          ht.pvalue(X1_hat, Y1_hat))
    print("Probability of X2 and Y1 under null hypothesis <=",
          ht.pvalue(X2_hat, Y1_hat))

# hypothesis testing with two halves of a real dataset
if 0:
    np.random.seed(0)
    
    # NOTE: these are local files which are unavailable on other machines
    
    # 3 cell line mixture
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    
    # pancreatic islets data
    #bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_bulk.csv"
    #scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_sc.csv"
    
    print("Reading full data matrices X and Y...")
    X = preprocessing.csvToMatrix(bulkfile)
    Y = preprocessing.csvToMatrix(scfile)
    
    N = X.shape[0]
    assert(Y.shape[0] == N)
    print("Number of genes:", N)
    print("Number of bulk samples:", X.shape[1])
    print("Number of single cells:", Y.shape[1])
    
    print("Shuffling X alongside Y...")
    indices = np.arange(N)
    np.random.shuffle(indices)
    X = X[indices]
    Y = Y[indices]
    
    print("Constructing half data matrices X1, X2, Y1, Y2...")
    halfway = int(N / 2)
    X1 = X[:halfway]
    X2 = X[halfway:]
    Y1 = Y[:halfway]
    Y2 = Y[halfway:]

    print("Are X1 and Y1 joint? Expect 0, receive",
          ht.residualsToCone(X1, Y1))
    print("Are X1 and Y2 joint? Expect > 0, receive",
          ht.residualsToCone(X1, Y2))
    print("Are X2 and Y1 joint? Expect > 0, receive",
          ht.residualsToCone(X2, Y1))
    print("Are X2 and Y2 joint? Expect 0, receive",
          ht.residualsToCone(X2, Y2))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          ht.pvalue(X1, Y1))
    print("Probability of X1 and Y2 under null hypothesis <=",
          ht.pvalue(X1, Y2))
    print("Probability of X2 and Y1 under null hypothesis <=",
          ht.pvalue(X2, Y1))
    print("Probability of X2 and Y2 under null hypothesis <=",
          ht.pvalue(X2, Y2))

# perform ZIFA on simulated data
if 0:
    np.set_printoptions(precision=3)
    
    N = 70
    M = 20
    L = 100
    K = 4
    lam = 0.1
    alpha = [ 1 for _ in range(K) ]
    #alpha = [9, 16, 0.3, 6]
    
    A = sims.randomA(N, K)
    X = sims.bulk(N, M, K, alpha, A=A)
    Y = sims.true_single_cell(N, L, K, alpha, A=A)
    print("noiseless bulk data:\n", X)
    print("noiseless single-cell data, no dropouts:\n", Y)
    sims.doubleExpDropouts(Y, lam)
    print("single-cell data after dropouts:\n", Y)
    
    Y = preprocessing.ZIFApreprocessing(Y)
    Z, params = ZIFA.fitModel(Y.T, K)
    print("ZIFA estimated latent positions:\n", Z.T)
    
# Plot description of variances of gene expressions in real scRNA-seq data
if 0:
    import preprocessing
    
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_islets_sc.csv"
    
    print("Reading single-cell data matrix Y...")
    Y = preprocessing.csvToMatrix(scfile)
    
    '''
    n = 400
    indices = set([ np.random.randint(Y.shape[0]) for _ in range(n) ])
    '''
    indices = np.arange(Y.shape[0])
    
    print("Computing variances for", len(indices), "random genes..." )
    V = [ np.var(Y[l]) for l in indices ]
    V.sort()
    
    print("Plotting variances...")
    plt.plot(np.arange(len(V)), V)
    plt.xlabel("gene")
    plt.ylabel("variance")

# assess solid angle of original and PCA cones for simulated single-cell 
# data
if 0:
    import simulations as sims
    from sklearn.decomposition import PCA
    
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X, Y = sims.simulateJointData()
    
    print("Performing PCA on single-cell data...")
    D = 206
    pca = PCA(n_components=D)
    Y_hat = pca.fit_transform(Y.T).T
    
    print("Scaled solid angle of original cone ~=", 
          util.scaledSolidAngle(Y))
    print("SSA of cone in PCA feature space ~=", 
          util.scaledSolidAngle(Y_hat))

# sample from 2d unit sphere first quadrant
if 0:
    points = util.uniformFromUnitSphere(2, 200)
    
    # reflect to positive orthant
    points = np.abs(points)
    
    # 2d case
    plt.scatter(
        points[:, 0],
        points[:, 1],
        )

# sample from 3d unit sphere positive orthant
if 0:
    points = util.uniformFromUnitSphere(3, 200)
    
    # reflect to positive orthant
    points = np.abs(points)
    
    # 3d case
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2]
        )
    
# Plot description of means and variances in real scRNA-seq data
if 0:
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    
    print("Reading single-cell data matrix Y...")
    Y = preprocessing.csvToMatrix(scfile)
    
    #n = 400
    #indices = set([ np.random.randint(Y.shape[0]) for _ in range(n) ])
    indices = np.arange(Y.shape[0])
    
    print("Computing expectation and variance for", len(indices), "genes...")
    E = np.array([ np.mean(Y[l]) for l in indices ]).astype(float)
    V = np.array([ np.var(Y[l]) for l in indices ]).astype(float)
    
    E = np.log(E + 1)
    V = np.log(V + 1)
    
    print("Plotting data...")
    
    fig, ax = plt.subplots()
    ax.set_xlim(left=0, right=max(E))
    ax.set_ylim(bottom=0, top=max(V))
    ax.scatter(E, V)
    ax.set_xlabel("log expectation")
    ax.set_ylabel("log variance")
    fig.tight_layout()
    plt.show()

# Plot description of mean and coefs of variation in real scRNA-seq data
if 0:
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    
    print("Reading single-cell data matrix Y...")
    Y = preprocessing.csvToMatrix(scfile)
    
    #n = 400
    #indices = set([ np.random.randint(Y.shape[0]) for _ in range(n) ])
    indices = np.arange(Y.shape[0])
    
    print("Computing expectation and var coef for", len(indices), "genes...")
    E = np.array([ np.mean(Y[l]) for l in indices ]).astype(float)
    V = np.array([ stats.variation(Y[l]) for l in indices ]).astype(float)
    
    E = np.log(E + 1)
    #V = np.log(V + 1)
    
    print("Plotting data...")
    
    fig, ax = plt.subplots()
    ax.set_xlim(left=0, right=max(E))
    ax.set_ylim(bottom=0, top=max(V))
    ax.scatter(E, V)
    ax.set_xlabel("log expectation")
    ax.set_ylabel("coef of variation")
    fig.tight_layout()
    plt.show()
    
# Plot description of mean and dropout rate in real scRNA-seq data
if 0:
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    
    print("Reading single-cell data matrix Y...")
    Y = preprocessing.csvToMatrix(scfile)
    
    #n = 400
    #indices = set([ np.random.randint(Y.shape[0]) for _ in range(n) ])
    indices = np.arange(Y.shape[0])
    
    print("Computing expectation and dropout rare for", 
          len(indices), "genes...")
    E = np.array([ np.mean(Y[n]) for n in indices ]).astype(float)
    V = np.array([ util.dropoutRate(Y[n]) for n in indices ]).astype(float)
    
    E = np.log(E + 1)
    
    print("Plotting data...")
    
    fig, ax = plt.subplots()
    ax.set_xlim(left=0, right=max(E))
    ax.set_ylim(bottom=0, top=max(V))
    ax.scatter(E, V)
    ax.set_xlabel("log expectation")
    ax.set_ylabel("dropout rate")
    fig.tight_layout()
    plt.show()

# assess solid angle of original and PCA cones for simulated single-cell 
# data
if 0:
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X, Y = sims.simulateJointData()
    
    print("Performing PCA on single-cell data...")
    D = 0.7
    pca = PCA(n_components=D)
    Y_hat = pca.fit_transform(Y.T).T
    
    print("Scaled solid angle of original cone ~=", 
          util.scaledSolidAngle(Y))
    print("SSA of cone in PCA feature space ~=", 
          util.scaledSolidAngle(Y_hat))

# sample from 2d unit sphere first quadrant
if 0:
    points = util.uniformFromUnitSphere(2, 200)
    
    # reflect to positive orthant
    points = np.abs(points)
    
    # 2d case
    plt.scatter(
        points[:, 0],
        points[:, 1],
        )

# sample from 3d unit sphere positive orthant
if 0:
    points = util.uniformFromUnitSphere(3, 200)
    
    # reflect to positive orthant
    points = np.abs(points)
    
    # 3d case
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2]
        )
    
# find background in simulated data... should be small
if 0:
    X, Y = sims.simulateJointData()
    bx, normx = util.inferBackground(X)
    by, normy = util.inferBackground(Y)

# assess similarity of profiles in islets and 3cl
if 0:
    three_cell = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3clwithisletssc.csv"
    islets = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_isletswith3clsc.csv"
    
    print("Loading datasets...")
    C = preprocessing.csvToMatrix(three_cell)
    I = preprocessing.csvToMatrix(islets) 
    Lc = C.shape[1]
    Li = C.shape[1]
    
    print("Normalizing columns...")
    C = normalize(C, norm='l1', axis=0)
    I = normalize(I, norm='l1', axis=0)
    
    print("Computing distances to mean...")
    c_mean = np.mean(C, axis=1)
    i_mean = np.mean(I, axis=1)
    
    print("Expectation of distance from mean in 3cl data =")
    c_devs = [ np.linalg.norm(C[:, l] - c_mean) for l in range(Lc) ]
    print(np.mean(c_devs))
    
    print("Expectation of distance from mean in islets data =")
    i_devs = [ np.linalg.norm(I[:, l] - i_mean) for l in range(Li) ]
    print(np.mean(i_devs))
    
    print("Distance between 3cl and islets means =", 
          np.linalg.norm(c_mean - i_mean))
    
    
    
    